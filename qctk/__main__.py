"""CLI entry point for qctk."""

import argparse
import sys

from logzero import logger

from qctk import __version__
from .config import StorageEngine
from .bam.extract import bam_extract_config_parser
from .fastq.kmers import fastq_kmers_config_parser
from .fastq.extract import fastq_extract_config_parser


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action="version", version="%(prog)s {}".format(__version__))

    parser.add_argument("--verbose", action="store_true", help="Enable verbose mode.")
    parser.add_argument("--quiet", action="store_true", help="Enable quiet mode.")

    parser.add_argument(
        "--reference", required=False, help="Path to reference FASTA (required for some commands)."
    )

    parser.add_argument("--storage-path", help="Path to storage.")
    parser.add_argument(
        "--storage-engine",
        default=StorageEngine.AUTO.value,
        choices=[e.value for e in StorageEngine],
        help="Path to storage.",
    )

    subparser = parser.add_subparsers(dest="command")
    bam_extract_config_parser(subparser)
    fastq_kmers_config_parser(subparser)
    fastq_extract_config_parser(subparser)

    args = parser.parse_args(argv)
    logger.info("Options: %s" % vars(args))
    if not args.command:
        logger.warn("No command provided!")
        parser.print_help()
        return 0
    else:
        return args.cmd(args)
