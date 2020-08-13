"""Models for storing VCF variant statistics information."""

import enum
import pathlib
import typing

import attr
import vcfpy


class Genotype(enum.Enum):
    #: Reference homozygous.
    REF = "0/0"
    #: Heterozygous alternative.
    HET = "0/1"
    #: Homozygous alternative.
    HOM = "1/1"


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


@attr.s(auto_attribs=True, frozen=True)
class VariantStats:
    """Statistics of a variant call.

    A missing genotype indicates a "no-call".  Missing coverage information indicates that the
    VCF file that this was generated from did not provide that information or a "no-call".
    """

    #: Genotype at the site.
    genotype: typing.Optional[Genotype]
    #: Total coverage at the site (ref + alt).
    total_cov: typing.Optional[int]
    #: Variant coverage at the site (alt).
    alt_cov: typing.Optional[int]


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
