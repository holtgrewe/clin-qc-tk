import os
import gzip

import pytest

from qctk.models import vcf, fastq


def test_kmer_info_headers():
    expected = ["genome_release", "chromosome", "position", "reference", "alternative", "ref_kmer"]
    assert fastq.KmerInfo.headers() == expected


def test_kmer_info_to_record(example_kmer_info_record, example_kmer_info):
    assert example_kmer_info.to_record() == example_kmer_info_record


def test_kmer_info_from_record(example_kmer_info_record, example_kmer_info):
    kmer_info = fastq.KmerInfo.from_record(example_kmer_info_record)
    assert kmer_info == example_kmer_info


def test_kmer_info_alt_kmer(example_kmer_info):
    assert example_kmer_info.alt_kmer == "NNNNNNNNNNTNNNNNNNNNN"


def test_kmer_info_read_kmer_infos_txt(fs, example_kmer_info_record, example_kmer_info):
    fs.create_file(
        "/test.txt",
        contents="\n".join(
            ["#" + "\t".join(fastq.KmerInfo.headers()), "\t".join(example_kmer_info_record),]
        )
        + "\n",
    )
    infos = fastq.read_kmer_infos(path="/test.txt")
    assert infos == [example_kmer_info]


def test_kmer_info_read_kmer_infos_gz(fs, example_kmer_info_record, example_kmer_info):
    fs.create_file(
        "/test.txt.gz",
        contents=gzip.compress(
            (
                "\n".join(
                    [
                        "#" + "\t".join(fastq.KmerInfo.headers()),
                        "\t".join(example_kmer_info_record),
                    ]
                )
                + "\n"
            ).encode("utf-8")
        ),
    )
    infos = fastq.read_kmer_infos(path="/test.txt.gz")
    assert infos == [example_kmer_info]


def test_kmer_info_write_kmer_infos_txt(fs, example_kmer_info):
    fastq.write_kmer_infos([example_kmer_info], path="/test.txt")
    assert os.path.exists("/test.txt")
    assert fastq.read_kmer_infos(path="/test.txt") == [example_kmer_info]


def test_kmer_info_write_kmer_infos_gz(fs, example_kmer_info):
    fastq.write_kmer_infos([example_kmer_info], path="/test.txt.gz")
    assert os.path.exists("/test.txt.gz")
    assert fastq.read_kmer_infos(path="/test.txt.gz") == [example_kmer_info]
