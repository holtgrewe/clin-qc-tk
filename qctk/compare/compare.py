"""Implementation of the "compare" command."""

import argparse
import itertools
import json
import math
import pathlib
import typing

import attr
import cattr
from logzero import logger
import numpy as np

from ..models import vcf
from ..config import DEFAULT_MIN_COV
from .config import CompareConfig


_TNumpySampleStats = typing.TypeVar("NumpySampleStats")


@attr.s(auto_attribs=True, frozen=True)
class NumpySampleStats:
    """Augmented SampleStats that has a numpy array."""

    #: Sample statistics information.
    stats: vcf.SampleStats
    #: Numpy array for fast comparison.
    arr: np.array

    @classmethod
    def from_sample_stats(
        cls, sample_stats: vcf.SampleStats, min_coverage: int
    ) -> _TNumpySampleStats:
        depths = []
        genotypes = []
        for site_stats in sample_stats.site_stats:
            variant_stats = site_stats.stats
            depths.append(variant_stats.total_cov or 0)
            genotypes.append(variant_stats.genotype or "0/0")
        return NumpySampleStats(
            stats=sample_stats,
            arr=np.array(
                [
                    [dp > min_coverage for dp in depths],
                    [gt != "0/0" for gt in genotypes],
                    [gt == "1/1" for gt in genotypes],
                ],
                dtype=bool,
            ),
        )


def _compute_stats(lhs: np.array, rhs: np.array) -> typing.Tuple[int, int, int, int, int, int]:
    """Compute values needed for relatedness."""
    # Obtain shortcuts...
    lhs_mask = lhs[0]
    lhs_is_alt = lhs[1]
    lhs_hom_alt = rhs[2]
    rhs_mask = rhs[0]
    rhs_is_alt = rhs[1]
    rhs_hom_alt = rhs[2]
    mask = lhs_mask & rhs_mask
    # Compute bit array
    lhs_ref = ~lhs_is_alt & mask
    rhs_ref = ~rhs_is_alt & mask
    arr_i = lhs_is_alt & ~lhs_hom_alt & mask
    arr_j = rhs_is_alt & ~rhs_hom_alt & mask
    arr_ij = arr_i & arr_j & mask
    # Compute partial terms
    het_i = np.count_nonzero(arr_i)
    het_j = np.count_nonzero(arr_j)
    het_i_j = np.count_nonzero(arr_ij)
    n_ibs_0 = np.count_nonzero(((lhs_ref & rhs_hom_alt) | (rhs_ref & lhs_hom_alt)) & mask)
    n_ibs_2 = np.count_nonzero(
        ((lhs_hom_alt & rhs_hom_alt) | (lhs_ref & rhs_ref) | (lhs_ref & rhs_ref)) & mask
    )
    n_ibs_1 = np.count_nonzero(mask) - n_ibs_0 - n_ibs_2
    # Return values.
    return (n_ibs_0, n_ibs_1, n_ibs_2, het_i, het_j, het_i_j)


def compute_similarity(left: NumpySampleStats, right: NumpySampleStats) -> vcf.SimilarityPair:
    assert left.stats.sample.name < right.stats.sample.name
    names = (left.stats.sample.name, right.stats.sample.name)
    vals = _compute_stats(left.arr, right.arr)
    return vcf.SimilarityPair(*itertools.chain(names, vals))


def compare_compare_run(config: CompareConfig) -> int:
    """Compare variant sites comparison.

    Entry point from configuration.
    """
    logger.info("Running compare")
    logger.info("Configuration: %s", config)

    if not config.common.storage_path:
        logger.error("You must provide --storage-path!")
        return 1
    storage_path = pathlib.Path(config.common.storage_path)
    if not storage_path.exists():
        logger.error("Storage path %s does not exist!", storage_path)
        return 1

    logger.info("Loading sample stats...")
    sample_stats = {}
    for path_stats in sorted(storage_path.glob("??/????/*-stats.json")):
        with path_stats.open("rt") as jsonf:
            stats = cattr.structure(json.load(jsonf), vcf.SampleStats)
            sample_stats[stats.sample.name] = NumpySampleStats.from_sample_stats(
                stats, config.min_cov
            )

    logger.info("Loading similarity stats...")
    similarity_pairs = {}
    for path_sim in sorted(storage_path.glob("??/????/*-sim.json")):
        with path_sim.open("rt") as jsonf:
            stats = cattr.structure(json.load(jsonf), typing.List[vcf.SimilarityPair])
            similarity_pairs[stats.key] = stats

    logger.info("Computing missing similarities...")
    for left in sample_stats.values():
        for right in sample_stats.values():
            key = (left.stats.sample.name, right.stats.sample.name)
            if key[0] < key[1] and key not in similarity_pairs:
                similarity_pairs[key] = compute_similarity(left, right)

    logger.info("Writing similarity updates...")
    out_sims = {}
    for pair in similarity_pairs.values():
        out_sims.setdefault(pair.sample_i, {})[pair.key] = pair
        out_sims.setdefault(pair.sample_j, {})[pair.key] = pair
    for name, sims in out_sims.items():
        sample_path = vcf.sample_path(config.common.storage_path, name)
        sims_path = pathlib.Path(str(sample_path).replace("-stats.json", "-sim.json"))
        with sims_path.open("wt") as jsonf:
            json.dump(cattr.unstructure(list(sorted(sims.values()))), jsonf)

    logger.info("All done. Have a nice day!")
    return 0


def compare_compare_main(args: argparse.Namespace) -> int:
    """Compare variant sites comparison.

    Entry point from argparse Namespace.
    """
    return compare_compare_run(CompareConfig.from_namespace(args))


def compare_compare_config_parser(subparsers: argparse._SubParsersAction) -> None:
    """Add command "compare" to argument parser."""
    parser = subparsers.add_parser("compare", help="Extract fingerprint information from BAM file.")
    parser.add_argument(
        "--hidden-cmd", dest="cmd", default=compare_compare_main, help=argparse.SUPPRESS
    )

    parser.add_argument(
        "--num-procs", type=int, default=1, help="Number of processors to use",
    )
    parser.add_argument(
        "--min-cov",
        type=int,
        default=DEFAULT_MIN_COV,
        help="Minimal depth to call a site covered, default: %d" % DEFAULT_MIN_COV,
    )
