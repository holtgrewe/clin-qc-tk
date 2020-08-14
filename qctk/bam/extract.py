import argparse
import os
import pathlib
import shlex
import subprocess
import tempfile
import typing

from logzero import logger
import pysam
import vcfpy

from .config import BamExtractConfig, DEFAULT_GENOME_RELEASE
from ..common import GenomeRelease, SITES_VCFS
from ..exceptions import SampleNameGuessingError
from ..models import vcf


#: Template for creating ``bcftools mpileup`` call.
TPL_PILEUP = r"bcftools mpileup -a AD,DP --threads 2 -I -R %(sites)s -f %(reference)s %(input_bam)s"

#: Template for creating ``bcftools call`` call.
TPL_CALL = r"bcftools call -c -Oz -o %(calls)s"


def calls_to_site_stats(
    config: BamExtractConfig, sites: typing.List[vcf.Site], calls_vcf: str
) -> typing.List[vcf.SiteStats]:
    """Merge calls with ``Site`` list to list of ``SiteStats``."""
    variant_stats = {site.short_notation: vcf.VariantStats() for site in sites}
    genome_release = sites[0].genome_release

    with vcfpy.Reader.from_path(calls_vcf) as reader:
        for record in reader:
            key = "-".join(map(str, (genome_release, record.CHROM, record.POS)))
            dp = record.INFO.get("DP")
            if record.ALT:
                ad = record.calls[0].data.get("AD")[1]
            else:
                ad = record.calls[0].data.get("AD")[0]
            gt = vcf.Genotype.from_value(record.calls[0].data.get("GT"))
            variant_stats[key] = vcf.VariantStats(genotype=gt, total_cov=dp, alt_cov=ad,)

    return [vcf.SiteStats(site=site, stats=variant_stats[site.short_notation],) for site in sites]


def call_sites(config: BamExtractConfig, path_sites_bed: str, path_bam: str, tmp_dir: str) -> str:
    """Perform variant calling at the given sites, write the results to a VCF file in ``tmp_dir``
    and return the path to the VCF file.
    """
    logger.info("Performing variant calling at sites")
    path_vcf = os.path.join(tmp_dir, "calls.vcf.gz")
    cmd_pileup = TPL_PILEUP % {
        "sites": path_sites_bed,
        "input_bam": path_bam,
        "reference": config.common.reference,
    }
    cmd_call = TPL_CALL % {"calls": path_vcf}
    logger.info("  mpileup: %s", " ".join(shlex.split(cmd_pileup)))
    logger.info("  call:    %s", " ".join(shlex.split(cmd_call)))
    p_pileup = subprocess.Popen(shlex.split(cmd_pileup), stdout=subprocess.PIPE)
    p_call = subprocess.Popen(shlex.split(cmd_call), stdin=p_pileup.stdout)
    p_call.wait()
    p_pileup.wait()
    return path_vcf


def write_sites_bed(
    config: BamExtractConfig,
    func_pref: typing.Callable[[str], str],
    sites: typing.List[vcf.Site],
    tmp_dir: str,
) -> str:
    """Write the given ``sites`` to a BED file inside ``tmp_dir`` and return its path."""
    path_bed = os.path.join(tmp_dir, "sites.bed")
    logger.info("Writing sites BED file to %s", path_bed)

    with open(path_bed, "wt") as bedf:
        for siteno, site in enumerate(sites):
            if not config.max_sites or site < config.max_sites:
                print(
                    "\t".join(map(str, (site.chromosome, site.position - 1, site.position,))),
                    file=bedf,
                )

    logger.info("Wrote %s sites", "{:,}".format(siteno))
    return path_bed


def identity(x: str) -> str:
    return x


def add_prefix(x: str) -> str:
    return "chr" + x


def strip_prefix(x: str) -> str:
    return x[3:]


def guess_sample_id(alif: pysam.AlignmentFile) -> str:
    result = None
    for line in alif.header.get("RG", []):
        if result:
            if line.get("SM", line.get("ID")) != result:  # pragma: no cover
                logger.error("Found more than one sample name in BAM file read group.")
                logger.error(
                    "Hint: check for problems with the data and override by specifying "
                    "--sample-id on the command line"
                )
                raise SampleNameGuessingError("Could not guess sample name")
        else:
            result = line.get("SM", line.get("ID"))

    if not result:  # pragma: no cover
        logger.error("Found no read group in BAM file.")
        logger.error(
            "Hint: check for problems with the data and override by specifying "
            "--sample-id on the command line"
        )
        raise SampleNameGuessingError("Could not guess sample name")
    else:
        return result


def _bam_extract_impl_for_file(
    config: BamExtractConfig, sites: typing.List[vcf.Site], path_bam: str
) -> None:
    """Perform site-wise extraction from one BAM file."""
    logger.info("Extracting from BAM file: %s", path_bam)
    sites_have_prefix = sites[0].chromosome.startswith("chr")
    with pysam.AlignmentFile(
        path_bam, mode="r", reference_filename=config.common.reference
    ) as alif:
        alis_have_prefix = alif.references[0].startswith("chr")
        if config.sample_id:
            sample_id = config.sample_id
        else:
            sample_id = guess_sample_id(alif)

    if sites_have_prefix and not alis_have_prefix:
        func_pref = strip_prefix
    elif not sites_have_prefix and alis_have_prefix:
        func_pref = add_prefix
    else:
        func_pref = identity

    with tempfile.TemporaryDirectory() as tmp_dir:
        path_sites_bed = write_sites_bed(config, func_pref, sites, str(tmp_dir))
        path_calls_vcf = call_sites(config, path_sites_bed, path_bam, str(tmp_dir))
        site_stats = calls_to_site_stats(config, sites, path_calls_vcf)
        sample_stats = vcf.SampleStats(sample=vcf.Sample(name=sample_id), site_stats=site_stats,)
        path_json = vcf.write_site_stats(sample_stats, config.common.storage_path, sample_id)
        logger.info("Wrote site-wise stats to %s", path_json)


def _bam_extract_impl(config: BamExtractConfig, sites: typing.List[vcf.Site]) -> None:
    """Perform site-wise extraction from BAM file(s)."""
    for path in config.input_files:
        _bam_extract_impl_for_file(config, sites, path)


def bam_extract_run(config: BamExtractConfig) -> int:
    """Extract site-wise information from BAM File using bcftools call.

    Entry point from configuration.
    """
    logger.info("Running bam-extract")
    logger.info("Configuration: %s", config)

    # TODO: storage plug and play
    pathlib.Path(config.common.storage_path).mkdir(parents=True, exist_ok=True)
    if not pathlib.Path(config.common.reference).exists():
        logger.error("The --reference path must exist!")
        return 1
    if not config.common.storage_path:
        logger.error("--storage-path must be provided!")
        return 1

    logger.info("Loading sites...")
    genome_release = (
        config.genome_release if config.genome_release else DEFAULT_GENOME_RELEASE.value
    )
    path_vcf = (
        config.sites_vcf
        if config.sites_vcf
        else SITES_VCFS[GenomeRelease.from_value(genome_release)]
    )
    sites = vcf.read_sites(path=path_vcf, genome_release=genome_release, max_sites=config.max_sites)

    logger.info("Analyzing BAM data...")
    _bam_extract_impl(config, sites)

    logger.info("All done. Have a nice day!")
    return 0


def bam_extract_main(args: argparse.Namespace) -> int:
    """Extract site-wise information from BAM File using bcftools call.

    Entry point from argparse Namespace.
    """
    return bam_extract_run(BamExtractConfig.from_namespace(args))


def bam_extract_config_parser(subparsers: argparse._SubParsersAction) -> None:
    """Add command "Bam-extract" to argument parser."""
    parser = subparsers.add_parser(
        "bam-extract", help="Extract fingerprint information from BAM file."
    )
    parser.add_argument(
        "--hidden-cmd", dest="cmd", default=bam_extract_main, help=argparse.SUPPRESS
    )

    parser.add_argument(
        "--sample-id",
        help="Optional sample identifier, overrides sample name read from read group in BAM file",
    )
    parser.add_argument(
        "--input-files",
        action="append",
        nargs="+",
        required=True,
        help="Path(s) to Bam files to process",
    )
    parser.add_argument(
        "--sites-vcf",
        help="Path to site kmers TSV file to use, if blank pick default by --genome-release",
    )
    parser.add_argument(
        "--genome-release",
        default=GenomeRelease.GRCH37.value,
        help="Name of the genome release to use (e.g., for built-in VCF file), default: %s"
        % GenomeRelease.GRCH37.value,
    )
