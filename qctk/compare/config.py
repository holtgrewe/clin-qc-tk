"""Configuration for the commands implemented in the ``compare`` module."""

import argparse
import attr
import cattr
import types
import typing

from ..config import CommonConfig, DEFAULT_MIN_COV

_TCompareConfig = typing.TypeVar("CompareConfig")


@attr.s(auto_attribs=True, frozen=True)
class CompareConfig:
    """Configuration for the ``compare`` command."""

    #: Common configuration.
    common: CommonConfig

    #: Degree of parallelism.
    num_procs: int = 1

    #: Minimal coverage.
    min_cov: int = DEFAULT_MIN_COV

    @classmethod
    def from_namespace(
        cls, ns: typing.Union[argparse.Namespace, types.SimpleNamespace]
    ) -> _TCompareConfig:
        return cattr.structure({"common": vars(ns), **vars(ns)}, cls)
