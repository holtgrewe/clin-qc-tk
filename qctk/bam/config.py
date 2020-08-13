"""Configuration for the commands implemented in the ``bam`` module."""

import attr
import cattr

from ..config import CommonConfig


class _BaseConfig:
    """Base class for the configuration."""

    #: Common configuration.
    common: CommonConfig

    def from_namespace(
        self, ns: typing.Union[argparse.Namespace, types.SimpleNamespace]
    ) -> _TCommonConfig:
        return cattr.structure({"common": vars(ns), **vars(ns)}, CommonConfig)


class BamStatsConfig(_BaseConfig):
    """Configuration for the ``bam-stats`` command."""


class BamExtractConfig(_BaseConfig):
    """Configuration for teh ``bam-extract`` command."""
