"""Configuration for the commands implemented in the ``fastq`` module."""

import argparse
import types
import typing

import attr
import cattr

from ..common import GenomeRelease, DEFAULT_GENOME_RELEASE, flatten_list
from ..config import CommonConfig, DEFAULT_KMER_LENGTH, DEFAULT_THRESHOLD


_TBaseConfig = typing.TypeVar("_BaseConfig")


@attr.s(auto_attribs=True, frozen=True)
class _BaseConfig:
    """Base class for the ``fastq-*`` configuration."""

    #: Common configuration.
    common: CommonConfig

    @classmethod
    def from_namespace(
        cls, ns: typing.Union[argparse.Namespace, types.SimpleNamespace]
    ) -> _TBaseConfig:
        return cattr.structure({"common": vars(ns), **vars(ns)}, cls)


@attr.s(auto_attribs=True, frozen=True)
class FastqKmersConfig(_BaseConfig):
    """Configuration for the ``fastq-kmers`` command."""

    #: The path to the output file with the k-mer information.
    output_tsv: str

    #: Maximal number of sites to read, no limit if ``None``.
    max_sites: typing.Optional[int] = None

    #: The k-mer length to extract for.
    kmer_length: int = DEFAULT_KMER_LENGTH

    #: The input path to the sites VCF file.  If ``None`` then the sites VCF file shipping with
    #: ``qctk`` will be used, selected by ``genome_release``.
    sites_vcf: typing.Optional[str] = None

    #: The genome release to select the sites VCF file for, by default GRCh37 will be used.
    genome_release: GenomeRelease = GenomeRelease.GRCH37


@attr.s(auto_attribs=True, frozen=True)
class FastqExtractConfig(_BaseConfig):
    """Configuration for the ``fastq-extract`` command."""

    #: ID of the sample under analysis.
    sample_id: str

    #: List with input FASTQ files to analyze.
    input_files: typing.List[str]

    #: The kmer info file to use for the analysis.  If ``None`` then the sites kmer file
    #: shipping with ``qctk`` will be used according to ``genome_release``.
    kmer_infos: typing.Optional[str]

    #: The genome release to select the sites VCF file for, by default GRCh37 will be used.
    genome_release: GenomeRelease = DEFAULT_GENOME_RELEASE

    #: The default threshold to use.
    threshold: float = DEFAULT_THRESHOLD

    @classmethod
    def from_namespace(
        cls, ns: typing.Union[argparse.Namespace, types.SimpleNamespace]
    ) -> _TBaseConfig:
        ns.input_files = flatten_list(ns.input_files)
        return cattr.structure({"common": vars(ns), **vars(ns)}, cls)
