"""Test for ``bam-extract``"""

import hashlib
import json

from qctk.config import CommonConfig, StorageEngine
from qctk.common import GenomeRelease
from qctk.bam import extract
from qctk.bam.config import BamExtractConfig
from qctk.__main__ import main


def test_bam_extract_run_smoke_test(tmp_path):
    path_storage = tmp_path / "storage"
    path_storage.mkdir()

    # Exercise the code.
    config = BamExtractConfig(
        common=CommonConfig(
            storage_path=str(path_storage), reference="tests/data/small/ref.fasta",
        ),
        input_files=["tests/data/small/reads_singleton1.bam",],
        sites_vcf="tests/data/small/vars_ontarget.vcf.gz",
        genome_release="test-small",
        sample_id="test-sample",
    )
    res = extract.bam_extract_run(config)

    sample_hash = hashlib.sha256("test-sample".encode("utf-8")).hexdigest()

    # Check results.
    assert res == 0
    output_path = path_storage / sample_hash[:2] / sample_hash[:4] / (sample_hash + "-stats.json")
    assert output_path.exists()
    with output_path.open("rt") as jsonf:
        data = json.load(jsonf)
    assert len(data) == 2
    assert len(data["site_stats"]) == 84
    assert data["site_stats"][1]["site"] == {
        "genome_release": "test-small",
        "chromosome": "contig",
        "position": 16570,
        "reference": "A",
        "alternative": "T",
    }
    assert data["site_stats"][1]["stats"] == {"genotype": "1/1", "total_cov": 76, "alt_cov": 73}


def test_bam_extract_via_args(mocker):
    mocker.patch.object(extract, "bam_extract_run")
    main(
        [
            "--storage-path",
            "/path/storage",
            "--reference",
            "/path/reference.fasta",
            "bam-extract",
            "--input-files",
            "/path/input.bam",
        ]
    )
    extract.bam_extract_run.assert_called_once_with(
        BamExtractConfig(
            common=CommonConfig(
                storage_path="/path/storage",
                verbose=False,
                quiet=False,
                storage_engine=StorageEngine.AUTO,
                reference="/path/reference.fasta",
            ),
            input_files=["/path/input.bam"],
            sites_vcf=None,
            sample_id=None,
            genome_release=GenomeRelease.GRCH37.value,
            max_sites=None,
        )
    )


def test_identity():
    assert extract.identity("chr1") == "chr1"
    assert extract.identity("1") == "1"


def test_add_prefix():
    assert extract.add_prefix("1") == "chr1"


def test_strip_prefix():
    assert extract.strip_prefix("chr1") == "1"
