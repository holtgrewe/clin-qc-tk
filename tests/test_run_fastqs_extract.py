"""Test for ``fastq-extract``"""

import hashlib
import json
import os
import sys

from qctk.config import CommonConfig
from qctk.fastq import extract
from qctk.fastq.config import FastqExtractConfig
from qctk.models import fastq


# TODO: run from command line / mock out fastq_extract_run


def test_fastq_extract_run_smoke_test(tmp_path, small_kmers_info_path):
    path_storage = tmp_path / "storage"
    path_storage.mkdir()

    # Exercise the code.
    config = FastqExtractConfig(
        common=CommonConfig(storage_path=path_storage, reference="tests/data/small/ref.fasta",),
        sample_id="test-sample",
        input_files=[
            "tests/data/small/reads_singleton1_1.fq.gz",
            "tests/data/small/reads_singleton1_2.fq.gz",
        ],
        kmer_infos=small_kmers_info_path,
        genome_release="test-small",
    )
    res = extract.fastq_extract_run(config)

    sample_hash = hashlib.sha256("test-sample".encode("utf-8")).hexdigest()

    # Check results.
    assert res == 0
    output_path = path_storage / sample_hash[:2] / (sample_hash + ".json")
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
    assert data[1]["stats"] == {"genotype": "1/1", "total_cov": 29, "alt_cov": 29}
