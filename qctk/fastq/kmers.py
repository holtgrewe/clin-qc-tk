"""Implementation of ``fastq-kmers``: build kmers TSV file from reference FASTA and a variants
VCF file.
"""

import argparse
import pathlib
import gzip
import typing

from logzero import logger
import pysam

from .config import FastqKmersConfig, DEFAULT_KMER_LENGTH, GenomeRelease, DEFAULT_GENOME_RELEASE
from ..models import vcf, fastq


#: Genome release to file name.
_SITE_FILES = {
    GenomeRelease.GRCH37: "sites.GRCh37.vcf.gz",
    GenomeRelease.GRCH38: "sites.hg38.vcf.gz",
}


def _fastq_kmers_run(
    config: FastqKmersConfig, sites: typing.Iterable[vcf.Site], outputf: typing.TextIO
) -> int:
    delta = config.kmer_length // 2
    with pysam.FastaFile(config.common.reference) as fastaf:
        print("#" + "\t".join(fastq.KmerInfo.headers()), file=outputf)
        for siteno, site in enumerate(sites):
            kmer = fastaf.fetch(site.chromosome, site.position - delta - 1, site.position + delta)
            if not kmer[delta] == site.reference:
                raise Exception(
                    "Expected %s:%d to be %s but is %s!"
                    % (site.chromosome, site.position, site.reference, kmer[delta])
                )
            arr = (
                config.genome_release,
                site.chromosome,
                site.position,
                site.reference,
                site.alternative,
                kmer,
            )
            print("\t".join(map(str, arr)), file=outputf)
    return 0


def fastq_kmers_run(config: FastqKmersConfig) -> int:
    """Extract k-mers from reference FASTA a sites list.

    Entry point from configuration.
    """
    logger.info("Running fastq-kmers")
    logger.info("Configuration: %s", config)

    if not pathlib.Path(config.output_tsv).parent.exists():
        logger.error("Path to output file %s does not exist!", config.output_tsv)
        return 1

    logger.info("Loading sites...")
    genome_release = (
        config.genome_release if config.genome_release else DEFAULT_GENOME_RELEASE.value
    )
    path_vcf = (
        config.sites_vcf
        if config.sites_vcf
        else _SITE_FILES[GenomeRelease.from_value(genome_release)]
    )
    sites = vcf.read_sites(path=path_vcf, genome_release=genome_release, max_sites=config.max_sites)

    logger.info("Generating kmers...")
    if config.output_tsv.endswith(".gz"):
        with gzip.open(config.output_tsv, "wt") as outputf:
            return _fastq_kmers_run(config, sites, outputf)
    else:
        with open(config.output_tsv, "wt") as outputf:
            return _fastq_kmers_run(config, sites, outputf)

    logger.info("All done. Have a nice day!")
    return 0


def fastq_kmers_main(args: argparse.Namespace) -> int:
    """Extract k-mers from reference FASTA a sites list.

    Entry point from argparse Namespace.
    """
    return fastq_kmers_run(FastqKmersConfig.from_namespace(args))


def fastq_kmers_config_parser(parser: argparse.ArgumentParser) -> None:
    """Add command "fastq-kmers" to argument parser."""
    parser.add_argument("--hidden-cmd", dest="cmd", default=fastq_kmers_run, help=argparse.SUPPRESS)

    parser.add_argument(
        "--output-tsv",
        required=True,
        help="Path to TSV file to write the k-mer infos to, .gz extension for compression",
    )
    parser.add_argument(
        "--kmer-length",
        default=DEFAULT_KMER_LENGTH,
        type=int,
        help="k-mer length to use, default: %d" % DEFAULT_KMER_LENGTH,
    )
    parser.add_argument(
        "--sites-vcf",
        help="Path to sites VCF file, by default the one shipping with qctk will be used",
    )
    parser.add_argument(
        "--genome-release",
        default=GenomeRelease.GRCH37.value,
        help="Name of the genome release to use (e.g., for built-in VCF file), default: %s"
        % GenomeRelease.GRCH37.value,
    )
    parser.add_argument(
        "--max-sites", type=int, help="Maximal number of sites to read, no limit if unset."
    )
