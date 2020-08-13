"""Configuration for the commands implemented in the ``fastq`` module."""

import argparse
import enum
import types
import typing

import attr
import cattr

from ..config import CommonConfig

#: Default k-mer length to use.
DEFAULT_KMER_LENGTH = 21

#: Default thershold to use.
DEFAULT_THRESHOLD = 0.1

_TGenomeRelease = typing.TypeVar("GenomeRelease")


class GenomeRelease(enum.Enum):
    #: GRCh37 release.
    GRCH37 = "GRCh37"
    #: GRCh38 release.
    GRCH38 = "GRCh38"

    @classmethod
    def from_value(cls, value: str) -> _TGenomeRelease:
        for release in cls:
            if release.value == value:
                return release
        raise ValueError("Could not get release for value %s" % value)


#: The default genome release.
DEFAULT_GENOME_RELEASE = GenomeRelease.GRCH37


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
    genome_release: GenomeRelease = GenomeRelease.GRCH37

    #: The default threshold to use.
    threshold: float = DEFAULT_THRESHOLD

    @classmethod
    def from_namespace(
        cls, ns: typing.Union[argparse.Namespace, types.SimpleNamespace]
    ) -> _TBaseConfig:
        ns.input_files = [item for sublist in ns.input_files for item in sublist]
        return cattr.structure({"common": vars(ns), **vars(ns)}, cls)
