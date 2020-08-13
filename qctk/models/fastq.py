"""Models for implementing the ``fastq-*`` commands."""

import attr
import gzip
import pathlib
import itertools
import typing

from .vcf import Site
from ..exceptions import KmerInfosIoError


_KmerInfo = typing.TypeVar("KmerInfo")


@attr.s(auto_attribs=True, frozen=True)
class KmerInfo:
    """Information about a k-mer that represets a marker site"""

    #: The corresponding site.
    site: Site
    #: The reference k-mer
    ref_kmer: str

    @property
    def alt_kmer(self) -> str:
        """Returns the alternative k-mer."""
        kmer_len = len(self.ref_kmer)
        return (
            self.ref_kmer[: kmer_len // 2]
            + self.site.alternative
            + self.ref_kmer[kmer_len // 2 + 1 :]
        )

    @staticmethod
    def headers() -> typing.List[str]:
        """Return header records."""
        return list(itertools.chain(map(lambda f: f.name, attr.fields(Site)), ["ref_kmer"]))

    def to_record(self) -> typing.List[str]:
        """Return to list of values, appropriate to write out to TSV file."""
        return list(map(str, itertools.chain(vars(self.site).values(), [self.ref_kmer])))

    @staticmethod
    def from_record(record: typing.List[str]) -> _KmerInfo:
        return KmerInfo(
            site=Site(
                genome_release=record[0],
                chromosome=record[1],
                position=int(record[2]),
                reference=record[3],
                alternative=record[4],
            ),
            ref_kmer=record[5],
        )


def read_kmer_infos(
    *,
    path: typing.Optional[typing.Union[str, pathlib.Path]] = None,
    stream: typing.Optional[typing.TextIO] = None
) -> typing.List[KmerInfo]:
    """Load kmer infos from the given file."""
    if bool(path) == bool(stream):
        raise ValueError("Exactly one of path and stream must be provided")  # pragma: no cover
    elif path:
        path_str = str(path)
        if path_str.endswith(".gz"):
            with gzip.open(path_str, "rt") as inputf:
                return read_kmer_infos(stream=inputf)
        else:
            with open(path_str, "rt") as inputf:
                return read_kmer_infos(stream=inputf)
    else:  # == elif path:
        result = []
        header = None
        for line in stream:
            arr = line.strip().split("\t")
            if not header:
                header = arr
                header[0] = header[0][1:]
                if header != KmerInfo.headers():
                    raise KmerInfosIoError(
                        "Invalid headers: %s, expected: %s" % (header, KmerInfo.headers())
                    )  # pragma: no cover
            else:
                result.append(KmerInfo.from_record(arr))
        return result


def write_kmer_infos(
    kmer_infos: typing.List[KmerInfo],
    *,
    path: typing.Optional[typing.Union[str, pathlib.Path]] = None,
    stream: typing.Optional[typing.TextIO] = None
):
    """Write kmer infos to the given file."""
    if bool(path) == bool(stream):
        raise ValueError("Exactly one of path and stream must be provided")  # pragma: no cover
    elif path:
        path_str = str(path)
        if path_str.endswith(".gz"):
            with gzip.open(path_str, "wt") as outputf:
                return write_kmer_infos(kmer_infos, stream=outputf)
        else:
            with open(path_str, "wt") as outputf:
                return write_kmer_infos(kmer_infos, stream=outputf)
    else:  # == elif path:
        print("#" + "\t".join(KmerInfo.headers()), file=stream)
        for kmer_info in kmer_infos:
            print("\t".join(kmer_info.to_record()), file=stream)
