"""Accession recognition and GEO FTP path construction."""

import pytest

from omixfetch.accession import Archive, geo_ftp_stub, geo_kind, parse_accession
from omixfetch.errors import AccessionError


@pytest.mark.parametrize(
    "raw,archive,normalised",
    [
        ("GSE164378", Archive.GEO, "GSE164378"),
        ("gse164378", Archive.GEO, "GSE164378"),
        ("  GSM5008737  ", Archive.GEO, "GSM5008737"),
        ("GPL24676", Archive.GEO, "GPL24676"),
    ],
)
def test_parse_accession(raw, archive, normalised):
    assert parse_accession(raw) == (archive, normalised)


@pytest.mark.parametrize(
    "raw",
    ["PRJNA672658", "prjeb1234", "PRJDB5", "SRP289408", "ERP1", "DRP2", "SRX9385726",
     "SRR12921178", "ERR1", "DRR2", "SAMN16565585", "SAMEA123", "SAMD9"],
)
def test_read_archive_accessions_route_to_sra(raw):
    assert parse_accession(raw) == (Archive.SRA, raw.upper())


@pytest.mark.parametrize("raw", ["", "   ", "PRJXX1", "SRZ123", "E-MTAB", "GSE", "12345"])
def test_parse_accession_rejects_unsupported(raw):
    with pytest.raises(AccessionError):
        parse_accession(raw)


@pytest.mark.parametrize("raw", ["E-MTAB-5061", "e-mtab-5061", "E-GEOD-12345"])
def test_arrayexpress_accessions_get_a_specific_message(raw):
    """ArrayExpress is recognised but deferred, so the error has to say so."""
    with pytest.raises(AccessionError, match="ArrayExpress"):
        parse_accession(raw)


@pytest.mark.parametrize(
    "accession,stub",
    [
        ("GSE164378", "GSE164nnn"),
        ("GSE1", "GSEnnn"),
        ("GSE999", "GSEnnn"),
        ("GSE1000", "GSE1nnn"),
        ("GSM5008737", "GSM5008nnn"),
        ("GPL24676", "GPL24nnn"),
    ],
)
def test_geo_ftp_stub(accession, stub):
    assert geo_ftp_stub(accession) == stub


def test_geo_kind():
    assert geo_kind("GSE1") == "series"
    assert geo_kind("GSM1") == "samples"
    assert geo_kind("GPL1") == "platforms"
