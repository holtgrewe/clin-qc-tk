"""Test for ``fastq-kmers``"""

from qctk.config import CommonConfig, StorageEngine, DEFAULT_THRESHOLD
from qctk.common import GenomeRelease
from qctk.fastq.config import FastqKmersConfig
from qctk.fastq import kmers
from qctk.models import fastq
from qctk.__main__ import main


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


def test_fastq_kmers_via_args(mocker):
    mocker.patch.object(kmers, "fastq_kmers_run")
    main(
        [
            "--reference",
            "/path/reference.fasta",
            "fastq-kmers",
            "--sites-vcf",
            "/path/sites.vcf.gz",
            "--output-tsv",
            "/path/output.tsv.gz",
        ]
    )
    kmers.fastq_kmers_run.assert_called_once_with(
        FastqKmersConfig(
            common=CommonConfig(
                storage_path=None,
                verbose=False,
                quiet=False,
                storage_engine=StorageEngine.AUTO,
                reference="/path/reference.fasta",
            ),
            sites_vcf="/path/sites.vcf.gz",
            output_tsv="/path/output.tsv.gz",
            max_sites=None,
            genome_release=GenomeRelease.GRCH37.value,
        )
    )
