"""Test for ``bam-extract``"""

import hashlib
import json

from qctk.config import CommonConfig
from qctk.bam import extract
from qctk.bam.config import BamExtractConfig


# TODO: run from command line / mock out bam_extract_run


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
    output_path = path_storage / sample_hash[:2] / sample_hash[:4] / (sample_hash + ".json")
    assert output_path.exists()
    with output_path.open("rt") as jsonf:
        data = json.load(jsonf)
    assert len(data) == 84
    assert data[1]["site"] == {
        "genome_release": "test-small",
        "chromosome": "contig",
        "position": 16570,
        "reference": "A",
        "alternative": "T",
    }
    assert data[1]["stats"] == {"genotype": "1/1", "total_cov": 76, "alt_cov": 73}
