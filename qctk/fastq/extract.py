"""Implementation of ``fastq-extract``: extract genotype and coverage information at marker
sites.
"""

import argparse
import hashlib
import itertools
import json
import pathlib
import typing

import cattr
from logzero import logger
import pysam

from .config import GenomeRelease, FastqExtractConfig, DEFAULT_GENOME_RELEASE, DEFAULT_THRESHOLD
from ..common import revcomp
from ..models.fastq import read_kmer_infos, KmerInfo
from ..models.vcf import SiteStats, VariantStats, Genotype, write_site_stats

#: Genome release to file name.
_KMER_FILES = {
    GenomeRelease.GRCH37: "kmers.GRCh37.tsv.gz",
    GenomeRelease.GRCH38: "kmers.hg38.tsv.gz",
}


class KmerCounter:
    """Helper class for iplementing k-mer counting."""

    def __init__(self, kmers):
        self.fwd = tuple(kmers)
        self.rev = tuple(map(revcomp, kmers))
        self.counts = {k: 0 for k in itertools.chain(self.fwd, self.rev)}

    def tally(self, kmer: str) -> None:
        if kmer in self.counts:
            self.counts[kmer] += 1


def _call_genotype(threshold: float, ref_depth: int, alt_depth: int) -> Genotype:
    total_depth = ref_depth + alt_depth
    if total_depth == 0:
        return None
    frac = alt_depth / total_depth
    if frac < threshold:
        return Genotype.REF
    elif frac > 1.0 - threshold:
        return Genotype.HOM
    else:
        return Genotype.HET


def _fastq_extract_impl(
    config: FastqExtractConfig, kmer_infos: typing.List[KmerInfo]
) -> pathlib.Path:
    kmers = []
    for kmer_info in kmer_infos:
        kmers.append(kmer_info.ref_kmer)
        kmers.append(kmer_info.alt_kmer)
    kmer_length = len(kmers[-1])

    counter = KmerCounter(kmers)
    for fastq_path in config.input_files:
        logger.debug("Processing FASTQ file: %s", fastq_path)
        with pysam.FastxFile(fastq_path) as inputf:
            for record in inputf:
                seq_length = len(record.sequence)
                for i in range(seq_length - kmer_length + 1):
                    kmer = record.sequence[i : i + kmer_length]
                    if len(kmer) != kmer_length:
                        raise ValueError(
                            "Kmer of invalid length (%d): %s, should be %d"
                            % (len(kmer), kmer, kmer_length)
                        )
                    else:
                        counter.tally(kmer)

    result = []
    for kmer_info in kmer_infos:
        ref_kmer = kmer_info.ref_kmer
        ref_depth = counter.counts[ref_kmer] + counter.counts[revcomp(ref_kmer)]
        alt_kmer = kmer_info.alt_kmer
        alt_depth = counter.counts[alt_kmer] + counter.counts[revcomp(alt_kmer)]
        result.append(
            SiteStats(
                site=kmer_info.site,
                stats=VariantStats(
                    genotype=_call_genotype(config.threshold, ref_depth, alt_depth),
                    total_cov=ref_depth + alt_depth,
                    alt_cov=alt_depth,
                ),
            )
        )

    return write_site_stats(result, config.common.storage_path, config.sample_id)


def fastq_extract_run(config: FastqExtractConfig) -> int:
    """Extract k-mers from reference FASTA a sites list.

    Entry point from configuration.
    """
    logger.info("Running fastq-extract")
    logger.info("Configuration: %s", config)

    # TODO: storage plug and play
    if not pathlib.Path(config.common.storage_path).exists():
        logger.error("Storage path does not exist: %s", config.common.storage_path)
        return 1

    logger.info("Loading kmers...")
    genome_release = (
        config.genome_release if config.genome_release else DEFAULT_GENOME_RELEASE.value
    )
    if config.kmer_infos:
        path_kmer_infos = config.kmer_infos
    else:
        path_kmer_infos = _KMER_FILES[GenomeRelease.from_value(genome_release)]
    kmer_infos = read_kmer_infos(path=path_kmer_infos)

    logger.info("Analyzing FASTQ data...")
    _fastq_extract_impl(config, kmer_infos)

    logger.info("All done. Have a nice day!")
    return 0


def fastq_extract_main(args: argparse.Namespace) -> int:
    """Extract k-mers from reference FASTA a sites list.

    Entry point from argparse Namespace.
    """
    return fastq_extract_run(FastqExtractConfig.from_namespace(args))


def fastq_extract_config_parser(parser: argparse.ArgumentParser) -> None:
    """Add command "fastq-extract" to argument parser."""
    parser.add_argument(
        "--hidden-cmd", dest="cmd", default=fastq_extract_run, help=argparse.SUPPRESS
    )

    parser.add_argument(
        "--sample-id", required=True, help="ID of the sample under analysis",
    )
    parser.add_argument(
        "--input-files",
        action="append",
        nargs="+",
        required=True,
        help="Path(s) to FASTQ files to process",
    )
    parser.add_argument(
        "--kmer-infos",
        help="Path to site kmers TSV file to use, if blank pick default by --genome-release",
    )
    parser.add_argument(
        "--genome-release",
        default=GenomeRelease.GRCH37.value,
        help="Name of the genome release to use (e.g., for built-in VCF file), default: %s"
        % GenomeRelease.GRCH37.value,
    )
    parser.add_argument(
        "--threshold",
        default=DEFAULT_THRESHOLD,
        type=float,
        help="Simple threshold to use for het/hom alt calls, default: %s" % DEFAULT_THRESHOLD,
    )
