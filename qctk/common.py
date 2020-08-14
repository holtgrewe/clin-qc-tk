import enum
import pathlib
import typing


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


#: Genome release to file name.
SITES_VCFS = {
    GenomeRelease.GRCH37: pathlib.Path(__file__).parent / "data" / "sites.GRCh37.vcf.gz",
    GenomeRelease.GRCH38: pathlib.Path(__file__).parent / "data" / "sites.hg38.vcf.gz",
}


#: Reverse complement lookup table.
_RC = {
    "c": "g",
    "C": "G",
    "g": "c",
    "G": "C",
    "t": "a",
    "T": "A",
    "a": "t",
    "A": "T",
}


def revcomp(nts: str) -> str:
    """Reverse-complement the given string.

    Characters not in "ACGTacgt" will be kept as identical.
    """
    return "".join(map(lambda x: _RC.get(x, x), nts))


def flatten_list(lst: typing.List[typing.List[typing.Any]]) -> typing.List[typing.Any]:
    """Flatten the given list of list of values into a list of values."""
    return [item for sublist in lst for item in sublist]
