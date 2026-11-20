"""File selection filters."""

import pytest

from omixfetch.filters import apply_filters, parse_size
from omixfetch.models import RemoteFile


@pytest.fixture
def files():
    return [
        RemoteFile(url="u1", filename="a_RAW.tar", size=10_000_000_000, kind="supplementary"),
        RemoteFile(url="u2", filename="b_counts.csv.gz", size=5_000, kind="supplementary"),
        RemoteFile(url="u3", filename="c.idf.txt", size=1_000, kind="metadata"),
        RemoteFile(url="u4", filename="d_matrix.mtx.gz", size=None, kind="supplementary"),
    ]


def test_include_glob(files):
    assert [f.filename for f in apply_filters(files, include=["*.gz"])] == [
        "b_counts.csv.gz",
        "d_matrix.mtx.gz",
    ]


def test_exclude_wins_over_include(files):
    result = apply_filters(files, include=["*.gz"], exclude=["*matrix*"])
    assert [f.filename for f in result] == ["b_counts.csv.gz"]


def test_kind_filter(files):
    assert [f.filename for f in apply_filters(files, kinds=["metadata"])] == ["c.idf.txt"]


def test_max_size_keeps_unknown_sizes(files):
    result = apply_filters(files, max_size=1_000_000)
    # The 10 GB tar is dropped; the unknown-size file is kept.
    assert [f.filename for f in result] == ["b_counts.csv.gz", "c.idf.txt", "d_matrix.mtx.gz"]


@pytest.mark.parametrize(
    "text,expected",
    [("500MB", 500 << 20), ("1.5GB", int(1.5 * (1 << 30))), ("2048", 2048), ("10K", 10 << 10)],
)
def test_parse_size(text, expected):
    assert parse_size(text) == expected
