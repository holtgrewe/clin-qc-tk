"""Configuration for the commands implemented in the ``vcf`` module."""

import attr
import cattr

from ..config import CommonConfig


class VcfExtractConfig(_BaseConfig):
    """Configuration for the ``vcf-extract`` command."""

    #: Common configuration.
    common: CommonConfig

    def from_namespace(
        self, ns: typing.Union[argparse.Namespace, types.SimpleNamespace]
    ) -> _TCommonConfig:
        return cattr.structure({"common": vars(ns), **vars(ns)}, CommonConfig)
