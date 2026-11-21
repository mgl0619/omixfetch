"""geofetch resolver backend (optional extra)."""

import pytest

from omixfetch.errors import AccessionError, ResolveError
from omixfetch.resolve import choose_resolver, resolve
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


def test_explicit_resolver_is_never_overridden():
    assert choose_resolver("GSE164378", "native") == "native"
    assert choose_resolver("GSM5008737", "geofetch") == "geofetch"


@pytest.mark.parametrize("accession", ["GSM5008737", "GPL24676"])
def test_auto_uses_native_for_non_series(accession):
    """geofetch's project API returns nothing for GSM and raises on GPL."""
    assert choose_resolver(accession, "auto") == "native"


def test_auto_follows_geofetch_availability(monkeypatch):
    import omixfetch.sources.geofetch_resolver as gr

    monkeypatch.setattr(gr, "geofetch_available", lambda: True)
    assert choose_resolver("GSE164378", "auto") == "geofetch"

    monkeypatch.setattr(gr, "geofetch_available", lambda: False)
    assert choose_resolver("GSE164378", "auto") == "native"


def test_resolved_by_is_recorded(monkeypatch):
    """A record says which code path produced it, for bug reports."""
    import sys

    # omixfetch/__init__.py rebinds the name "resolve" to the function, so the
    # module itself has to come from sys.modules rather than attribute access.
    rv = sys.modules["omixfetch.resolve"]

    monkeypatch.setattr(rv, "choose_resolver", lambda acc, req="auto": "native")
    captured = {}

    def fake(accession, session, **kwargs):
        from omixfetch.models import Dataset, DatasetMetadata

        captured["kwargs"] = kwargs
        return Dataset(metadata=DatasetMetadata(accession=accession, archive="geo"))

    monkeypatch.setitem(rv._RESOLVERS, rv.Archive.GEO, fake)
    dataset = rv.resolve("GSE1")
    assert dataset.metadata.resolved_by == "native"
