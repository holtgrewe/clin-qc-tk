"""Configuration for the commands implemented in the ``bam`` module."""

import argparse
import attr
import cattr
import types
import typing

from ..common import GenomeRelease, DEFAULT_GENOME_RELEASE, flatten_list
from ..config import CommonConfig

_TBaseConfig = typing.TypeVar("_BaseConfig")


@attr.s(auto_attribs=True, frozen=True)
class _BaseConfig:
    """Base class for the configuration."""

    #: Common configuration.
    common: CommonConfig

    @classmethod
    def from_namespace(
        cls, ns: typing.Union[argparse.Namespace, types.SimpleNamespace]
    ) -> _TBaseConfig:
        ns.input_files = flatten_list(ns.input_files)
        return cattr.structure({"common": vars(ns), **vars(ns)}, cls)


@attr.s(auto_attribs=True, frozen=True)
class BamStatsConfig(_BaseConfig):
    """Configuration for the ``bam-stats`` command."""


@attr.s(auto_attribs=True, frozen=True)
class BamExtractConfig(_BaseConfig):
    """Configuration for teh ``bam-extract`` command."""

    #: Input BAM file to analyze.
    input_files: typing.List[str]

    #: The sites VCF file to use for the analysis.  If ``None`` then the sites kmer file
    #: shipping with ``qctk`` will be used according to ``genome_release``.
    sites_vcf: typing.Optional[str] = None

    #: Optional sample ID, overrides the one from the BAM file's read group.
    sample_id: typing.Optional[str] = None

    #: The genome release to select the sites VCF file for, by default GRCh37 will be used.
    genome_release: str = DEFAULT_GENOME_RELEASE.value

    #: The maximal number of sites to process, ``None``/``0`` for no limit.
    max_sites: typing.Optional[int] = None
