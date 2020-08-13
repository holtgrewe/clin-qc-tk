"""Test for ``fastq-kmers``"""

import os

from qctk.fastq.config import FastqKmersConfig
from qctk.config import CommonConfig
from qctk.fastq import kmers
from qctk.models import fastq


# TODO: run from command line / mock out fastq_kmers_run


def test_fastq_kmers_run_smoke_test(tmp_path):
    path_out = tmp_path / "out" / "out.tsv"
    path_out.parent.mkdir()

    # Exercise the code.
    config = FastqKmersConfig(
        common=CommonConfig(
            storage_path=tmp_path / "storage", reference="tests/data/small/ref.fasta",
        ),
        output_tsv=str(path_out),
        sites_vcf="tests/data/small/vars_ontarget.vcf.gz",
        genome_release="test-small",
    )
    res = kmers.fastq_kmers_run(config)

    # Check results.
    assert res == 0
    assert path_out.exists()
    lines = path_out.read_text().split("\n")
    assert len(lines) == 86
    assert lines[0] == "#" + "\t".join(fastq.KmerInfo.headers())
    assert lines[1] == "test-small\tcontig\t16556\tT\tA\tATTTTCCTCATTGAAGATATT"
    assert lines[-2] == "test-small\tcontig\t82577\tC\tG\tCCCGCCATACCGGCCGGGGCG"
    assert lines[-1] == ""
