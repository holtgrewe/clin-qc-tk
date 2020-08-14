"""Test for ``fastq-extract``"""

import hashlib
import json

from qctk.config import CommonConfig, StorageEngine, DEFAULT_THRESHOLD
from qctk.common import GenomeRelease
from qctk.fastq import extract
from qctk.fastq.config import FastqExtractConfig
from qctk.__main__ import main


# TODO: run from command line / mock out fastq_extract_run


def test_fastq_extract_run_smoke_test(tmp_path, small_kmers_info_path):
    path_storage = tmp_path / "storage"
    path_storage.mkdir()

    # Exercise the code.
    config = FastqExtractConfig(
        common=CommonConfig(
            storage_path=str(path_storage), reference="tests/data/small/ref.fasta",
        ),
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
    assert data[1]["stats"] == {"genotype": "1/1", "total_cov": 29, "alt_cov": 29}


def test_fastq_extract_via_args(mocker):
    mocker.patch.object(extract, "fastq_extract_run")
    main(
        [
            "--storage-path",
            "/path/storage",
            "fastq-extract",
            "--sample-id",
            "test-sample",
            "--input-files",
            "/path/input.fastq.gz",
        ]
    )
    extract.fastq_extract_run.assert_called_once_with(
        FastqExtractConfig(
            common=CommonConfig(
                storage_path="/path/storage",
                verbose=False,
                quiet=False,
                storage_engine=StorageEngine.AUTO,
                reference=None,
            ),
            input_files=["/path/input.fastq.gz"],
            sample_id="test-sample",
            kmer_infos=None,
            genome_release=GenomeRelease.GRCH37.value,
            threshold=DEFAULT_THRESHOLD,
        )
    )
