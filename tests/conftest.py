import pytest

from qctk.config import CommonConfig
from qctk.fastq import kmers
from qctk.fastq.config import FastqKmersConfig
from qctk.models import vcf, fastq


@pytest.fixture
def example_site():
    return vcf.Site(
        genome_release="GRCh37", chromosome="1", position=123, reference="C", alternative="T",
    )


@pytest.fixture
def example_kmer_info(example_site):
    return fastq.KmerInfo(site=example_site, ref_kmer="NNNNNNNNNNANNNNNNNNNN")


@pytest.fixture
def example_kmer_info_record():
    return [
        "GRCh37",
        "1",
        "123",
        "C",
        "T",
        "NNNNNNNNNNANNNNNNNNNN",
    ]


@pytest.fixture
def small_kmers_info_path(tmp_path):
    path_out = tmp_path / "kmers" / "kmers.tsv"
    path_out.parent.mkdir()

    config = FastqKmersConfig(
        common=CommonConfig(
            storage_path=tmp_path / "storage", reference="tests/data/small/ref.fasta",
        ),
        output_tsv=str(path_out),
        sites_vcf="tests/data/small/vars_ontarget.vcf.gz",
        genome_release="test-small",
    )
    res = kmers.fastq_kmers_run(config)

    assert res == 0
    return str(path_out)
