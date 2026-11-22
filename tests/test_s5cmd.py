"""s5cmd backend: S3 URL translation and batch command construction."""

import pytest

from omixfetch.backends.s5cmd import (
    S5cmdBackend,
    render_run_file,
    s3_uri,
    s5cmd_available,
)
from omixfetch.models import RemoteFile


@pytest.mark.parametrize(
    "url,expected",
    [
        (
            "https://sra-pub-src-2.s3.amazonaws.com/SRR16036138/x.bam.1",
            "s3://sra-pub-src-2/SRR16036138/x.bam.1",
        ),
        (
            "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR390728/SRR390728",
            "s3://sra-pub-run-odp/sra/SRR390728/SRR390728",
        ),
        (
            "https://sra-pub-src-1.s3.us-east-1.amazonaws.com/a/b.bam",
            "s3://sra-pub-src-1/a/b.bam",
        ),
        (
            "https://s3.us-east-1.amazonaws.com/bucket/key/file.bam",
            "s3://bucket/key/file.bam",
        ),
    ],
)
def test_s3_uri_translates_both_url_styles(url, expected):
    assert s3_uri(url) == expected


@pytest.mark.parametrize(
    "url",
    [
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/SRR1039514_1.fastq.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1/suppl/x.gz",
        "",
    ],
)
def test_s3_uri_returns_none_for_non_s3(url):
    assert s3_uri(url) is None


def test_render_run_file_uses_layout_targets(tmp_path):
    f = RemoteFile(
        url="https://sra-pub-src-2.s3.amazonaws.com/SRR1/x.bam",
        filename="x.bam",
    )
    f.relpath = "GSM1/x.bam"
    text = render_run_file([f], tmp_path, concurrency=8)
    assert text.startswith("cp --concurrency 8 's3://sra-pub-src-2/SRR1/x.bam'")
    assert str(tmp_path / "GSM1" / "x.bam") in text


def test_render_run_file_skips_non_s3_entries(tmp_path):
    files = [
        RemoteFile(url="https://sra-pub-src-2.s3.amazonaws.com/a/x.bam", filename="x.bam"),
        RemoteFile(url="https://ftp.sra.ebi.ac.uk/vol1/y.gz", filename="y.gz"),
    ]
    assert len(render_run_file(files, tmp_path).strip().splitlines()) == 1


def test_render_run_file_empty(tmp_path):
    assert render_run_file([], tmp_path) == ""


@pytest.mark.skipif(not s5cmd_available(), reason="s5cmd not installed")
def test_build_command_carries_parallelism(tmp_path):
    cmd = S5cmdBackend(concurrent=8, connections=16).build_command(tmp_path / "f.txt")
    assert "--no-sign-request" in cmd
    assert "--numworkers" in cmd and "8" in cmd
    assert cmd[-2] == "run"
