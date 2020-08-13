"""Common code for configuration in ``qctk``."""

import enum
import types
import typing

import argparse
import attr
import cattr


#: Enumeration of the supported storage engines.
class StorageEngine(enum.Enum):
    #: Auto-detect storage engine by path's extension.
    AUTO = "auto"
    #: File system storage engine (storage path is a base directory).
    FS = "fs"


#: Forward declration of ``CommonConfig``.
_TCommonConfig = typing.TypeVar("CommonConfig")


@attr.s(auto_attribs=True, frozen=True)
class CommonConfig:
    """Configuration for all ``qctk`` commands."""

    #: Path to the storage directory/file.
    storage_path: typing.Optional[str]

    #: Enables the "verbose" mode (takes precedence over ``quiet``).
    verbose: bool = False

    #: Enables the "quiet" mode.
    quiet: bool = False

    #: The storage engine to use.
    storage_engine: StorageEngine = StorageEngine.AUTO

    #: Path to the FAI-indexed reference FASTA file.  This is required when using CRAM files
    #: but also for a couple of other commands.
    reference: typing.Optional[str] = None

    @property
    def is_verbose(self) -> bool:
        """Returns whether verbose mode is enabled, for symmatrie with ``is_quiet``."""
        return self.verbose

    @property
    def is_quiet(self) -> bool:
        """Returns whether quiet mode is enabled (verbose takes precedence)."""
        return not self.verbose and self.quiet

    def from_namespace(
        self, ns: typing.Union[argparse.Namespace, types.SimpleNamespace]
    ) -> _TCommonConfig:
        return cattr.structure(vars(ns), CommonConfig)
