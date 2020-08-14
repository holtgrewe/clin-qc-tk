"""Exception classes for ``qctk``."""


class QctkException(Exception):
    """Base exception class."""


class KmerInfosIoError(Exception):
    """Raised when kmer infos could not be loaded."""


class SampleNameGuessingError(Exception):
    """Raised when sample name guessing failed."""
