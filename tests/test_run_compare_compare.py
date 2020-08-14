"""Test for ``compare``"""

import pathlib
import json

import cattr

from qctk.config import CommonConfig, StorageEngine
from qctk.compare.config import CompareConfig
from qctk.compare import compare
from qctk.models import vcf
from qctk.__main__ import main


def test_compare_compare_run_smoke_test(tmp_path):
    sample_paths = []
    sample_ids = [
        "singleton1",
        "singleton2",
    ]
    for i, sample_id in enumerate(sample_ids):
        sample_path = vcf.sample_path(str(tmp_path / "storage"), sample_id)
        sample_paths.append(sample_path)
        sample_path.parent.mkdir(exist_ok=True, parents=True)
        sample_stats = vcf.SampleStats(
            sample=vcf.Sample(name=sample_id),
            site_stats=[
                vcf.SiteStats(
                    site=vcf.Site(
                        genome_release="GRCh37",
                        chromosome="1",
                        position=j + 1,
                        reference="A",
                        alternative="T",
                    ),
                    stats=vcf.VariantStats(
                        genotype=vcf.Genotype.HET if i == j else vcf.Genotype.REF,
                        total_cov=30,
                        alt_cov=15 if i == j else 0,
                    ),
                )
                for j in range(2)
            ],
        )
        with sample_path.open("wt") as samplef:
            json.dump(cattr.unstructure(sample_stats), samplef)

    # Exercise the code.
    config = CompareConfig(
        common=CommonConfig(
            storage_path=tmp_path / "storage", reference="tests/data/small/ref.fasta",
        ),
    )
    res = compare.compare_compare_run(config)

    # Check results.
    assert res == 0
    for sample_path in sample_paths:
        p = pathlib.Path(str(sample_path).replace("-stats.json", "-sim.json"))
        assert p.exists()
        with p.open("rt") as jsonf:
            vals = json.load(jsonf)
            assert vals == [
                {
                    "sample_i": "singleton1",
                    "sample_j": "singleton2",
                    "n_ibs0": 0,
                    "n_ibs1": 2,
                    "n_ibs2": 0,
                    "het_i": 2,
                    "het_j": 2,
                    "het_i_j": 2,
                }
            ]


def test_compare_compare_via_args(mocker):
    mocker.patch.object(compare, "compare_compare_run")
    main(
        ["--storage-path", "/path/storage", "compare",]
    )
    compare.compare_compare_run.assert_called_once_with(
        CompareConfig(
            common=CommonConfig(
                storage_path="/path/storage",
                verbose=False,
                quiet=False,
                storage_engine=StorageEngine.AUTO,
                reference=None,
            ),
            num_procs=1,
        )
    )
