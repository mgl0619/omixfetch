"""SOFT parsing and GEO record construction."""

import requests

from omixfetch.http import get_text
from omixfetch.sources.geo import (
    _series_metadata,
    _sample_metadata,
    _supplementary_files,
    parse_soft,
    series_matrix_url,
)


def test_parse_soft_splits_entities(series_soft, samples_soft):
    series = parse_soft(series_soft)
    assert len(series) == 1
    entity, accession, fields = series[0]
    assert (entity, accession) == ("SERIES", "GSE000001")
    assert fields["Series_title"] == ["A small test series"]
    # Repeated keys are preserved rather than overwritten.
    assert len(fields["Series_summary"]) == 2
    assert len(fields["Series_sample_id"]) == 2

    samples = parse_soft(samples_soft)
    assert [a for _, a, _ in samples] == ["GSM000001", "GSM000002"]


def test_series_metadata_fields(series_soft):
    _, _, fields = parse_soft(series_soft)[0]
    m = _series_metadata("GSE000001", fields)
    assert m.title == "A small test series"
    assert m.summary == "First summary paragraph. Second summary paragraph."
    assert m.n_samples == 2
    assert m.organisms == ["Homo sapiens"]
    assert m.pubmed_ids == ["33456401"]
    assert m.release_date == "Jan 15 2021"
    assert m.relations["BioProject"].endswith("PRJNA000001")
    assert m.contact["institute"] == "Test Institute"
    assert len(m.contributors) == 2


def test_sample_characteristics_become_key_values(samples_soft):
    records = parse_soft(samples_soft)
    second = _sample_metadata("GSM000002", records[1][2])
    assert second.characteristics == {
        "cell type": "PBMC",
        "treatment": "stimulated",
        "timepoint": "24h",
    }
    assert second.library_strategy == "RNA-Seq"


def test_supplementary_urls_are_https_and_skip_none(series_soft, samples_soft):
    _, _, series_fields = parse_soft(series_soft)[0]
    files = _supplementary_files(series_fields, "Series", "GSE000001")
    assert len(files) == 2
    assert all(f.url.startswith("https://ftp.ncbi.nlm.nih.gov/") for f in files)
    assert {f.filename for f in files} == {"GSE000001_RAW.tar", "GSE000001_counts.csv.gz"}

    # "NONE" placeholders must not become downloads.
    _, _, sample_fields = parse_soft(samples_soft)[1]
    sample_files = _supplementary_files(sample_fields, "Sample", "GSM000002")
    assert len(sample_files) == 1


def test_series_matrix_url():
    assert series_matrix_url("GSE164378") == (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164378/matrix/"
        "GSE164378_series_matrix.txt.gz"
    )


def test_get_text_forces_utf8_when_charset_is_absent(monkeypatch):
    """NCBI serves UTF-8 without a charset; requests would guess Latin-1."""

    class FakeResponse:
        ok = True
        status_code = 200
        headers = {"Content-Type": "text/plain"}
        encoding = "ISO-8859-1"
        _raw = "!Series_summary = \u201cweighted-nearest neighbor\u201d".encode("utf-8")

        @property
        def text(self):
            return self._raw.decode(self.encoding)

    session = requests.Session()
    monkeypatch.setattr(session, "get", lambda *a, **k: FakeResponse())
    assert "\u201cweighted-nearest neighbor\u201d" in get_text(session, "https://example.org")
