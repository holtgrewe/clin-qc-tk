"""Models for storing VCF variant statistics information."""

import enum
import hashlib
import pathlib
import typing

import attr
import cattr
import json
from logzero import logger
import vcfpy


_TGenotype = typing.TypeVar("Genotype")


class Genotype(enum.Enum):
    #: Reference homozygous.
    REF = "0/0"
    #: Heterozygous alternative.
    HET = "0/1"
    #: Homozygous alternative.
    HOM = "1/1"

    @classmethod
    def from_value(cls, value: str) -> _TGenotype:
        val = value.replace("|", "/")
        for release in cls:
            if release.value == val:
                return release
        raise ValueError("Could not get release for value %s" % value)


@attr.s(auto_attribs=True, frozen=True)
class Site:
    """Define a single site."""

    #: The genome release
    genome_release: str
    #: The contig/chromosome name.
    chromosome: str
    #: 1-based position of the site.
    position: int
    #: Reference base string in VCF notation.
    reference: str
    #: Alternative base string in VCF notation.
    alternative: str

    def with_prefix(self, prefix):
        """Return ``Site`` having the ``"chr"`` prefix or not."""
        if prefix and not self.chromosome.startswith("chr"):
            return attr.evolve(self, chromosome="chr" + self.chromosome)
        elif not prefix and self.chromosome.startswith("chr"):
            return attr.evolve(self, chromosome=self.chromosome[3:])

    @property
    def short_notation(self) -> str:
        return "-".join(map(str, [self.genome_release, self.chromosome, self.position]))


@attr.s(auto_attribs=True, frozen=True)
class VariantStats:
    """Statistics of a variant call.

    A missing genotype indicates a "no-call".  Missing coverage information indicates that the
    VCF file that this was generated from did not provide that information or a "no-call".
    """

    #: Genotype at the site.
    genotype: typing.Optional[Genotype] = None
    #: Total coverage at the site (ref + alt).
    total_cov: typing.Optional[int] = None
    #: Variant coverage at the site (alt).
    alt_cov: typing.Optional[int] = None


@attr.s(auto_attribs=True, frozen=True)
class SiteStats:
    """Variant statistics at a site."""

    #: Site
    site: Site
    #: Variant statistics.
    stats: VariantStats


def read_sites(
    *,
    path: typing.Optional[typing.Union[str, pathlib.Path]] = None,
    stream: typing.Optional[vcfpy.Reader] = None,
    genome_release: typing.Optional[str] = None,
    max_sites: typing.Optional[int] = None,
) -> typing.List[Site]:
    """Load sites from the given VCF file."""
    if not genome_release:
        raise ValueError("genome_release must be given")  # pragma: no cover
    if bool(path) == bool(stream):
        raise ValueError("Exactly one of path and stream must be provided")  # pragma: no cover
    else:
        if path:
            reader = vcfpy.Reader.from_path(path)
        else:
            reader = vcfpy.Reader.from_stream(stream)

        result = []
        for record in reader:
            for lineno, record in enumerate(reader):
                if not max_sites or lineno < max_sites:
                    result.append(
                        Site(
                            genome_release=genome_release,
                            chromosome=record.CHROM,
                            position=record.POS,
                            reference=record.REF,
                            alternative=record.ALT[0].value,
                        )
                    )
                else:
                    break

    return result


def write_site_stats(site_stats: typing.List[SiteStats], storage_path: str, sample_id: str) -> str:
    """Write site stats to the storage path and return path to JSON."""
    # TODO: pluggable storage
    sample_hash = hashlib.sha256(sample_id.encode("utf-8")).hexdigest()
    output_path = (
        pathlib.Path(storage_path) / sample_hash[:2] / sample_hash[:4] / (sample_hash + ".json")
    )
    output_path.parent.mkdir(parents=True)
    logger.info("Writing results to %s", output_path)
    with output_path.open("wt") as outputf:
        json.dump(cattr.unstructure(site_stats), outputf)
    return output_path
