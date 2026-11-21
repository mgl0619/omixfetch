"""geofetch resolver backend (optional extra)."""

import pytest

from omixfetch.errors import AccessionError, ResolveError
from omixfetch.resolve import resolve
from omixfetch.sources.geofetch_resolver import (
    _cell,
    _https,
    _size,
    geofetch_available,
)


@pytest.mark.parametrize(
    "raw,expected",
    [
        (
            "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1/suppl/x.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1/suppl/x.gz",
        ),
        ("https://example.org/x.gz", "https://example.org/x.gz"),
        ("", ""),
    ],
)
def test_https_rewrites_geo_ftp_urls(raw, expected):
    assert _https(raw) == expected


def test_cell_maps_pandas_missing_values_to_empty_string():
    assert _cell({"a": "nan"}, "a") == ""
    assert _cell({"a": None}, "a") == ""
    assert _cell({"a": "  x  "}, "a") == "x"
    assert _cell({}, "missing") == ""


def test_size_parses_floats_and_blanks():
    assert _size({"file_size": "683095"}) == 683095
    assert _size({"file_size": "683095.0"}) == 683095
    assert _size({"file_size": ""}) is None
    assert _size({"file_size": "huge"}) is None


def test_unknown_resolver_is_rejected():
    with pytest.raises(AccessionError, match="Unknown resolver"):
        resolve("GSE164378", resolver="nonsense")


@pytest.mark.skipif(geofetch_available(), reason="geofetch is installed")
def test_missing_geofetch_gives_an_actionable_error():
    with pytest.raises(ResolveError, match="omixfetch\\[geofetch\\]"):
        resolve("GSE164378", resolver="geofetch")


@pytest.mark.skipif(not geofetch_available(), reason="geofetch not installed")
def test_geofetch_backend_rejects_non_series_accessions():
    with pytest.raises(ResolveError, match="GSE series only"):
        resolve("GSM5008737", resolver="geofetch")
