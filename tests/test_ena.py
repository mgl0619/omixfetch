"""ENA raw-read resolution."""

import pytest

from omixfetch.errors import ResolveError
from omixfetch.models import Dataset, DatasetMetadata, SampleMetadata
from omixfetch.sources.ena import (
    attach_reads,
    parse_filereport,
    project_accessions,
    _as_url,
)

FILEREPORT = (
    "run_accession\texperiment_accession\tsample_accession\tsecondary_sample_accession\t"
    "study_accession\tlibrary_strategy\tlibrary_source\tlibrary_selection\t"
    "instrument_model\tread_count\tbase_count\tfastq_ftp\tfastq_bytes\tfastq_md5\t"
    "submitted_ftp\tsubmitted_bytes\tsubmitted_md5\n"
    "SRR1039514\tSRX384345\tSAMN02422669\tSRS508568\tPRJNA229998\tRNA-Seq\t"
    "TRANSCRIPTOMIC\tcDNA\tIllumina HiSeq 2000\t40021727\t4042594427\t"
    "ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/004/SRR1039514/SRR1039514_1.fastq.gz;"
    "ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/004/SRR1039514/SRR1039514_2.fastq.gz\t"
    "2377036173;2376896820\tae35fe0ce13badacc48c65717e811528;9ef4fe59d6378c513f933e24d12f6047\t"
    "\t\t\n"
)


def test_project_accessions_from_geo_relation_urls():
    relations = {
        "BioProject": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA229998",
        "SRA": "https://www.ncbi.nlm.nih.gov/sra?term=SRP033351",
    }
    # BioProject is tried first because GEO records it more consistently.
    assert project_accessions(relations) == ["PRJNA229998", "SRP033351"]


def test_project_accessions_empty_when_no_links():
    assert project_accessions({"Other": "https://example.org/nothing"}) == []


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("ftp.sra.ebi.ac.uk/vol1/x.gz", "https://ftp.sra.ebi.ac.uk/vol1/x.gz"),
        ("ftp://ftp.sra.ebi.ac.uk/vol1/x.gz", "https://ftp.sra.ebi.ac.uk/vol1/x.gz"),
        ("https://ftp.sra.ebi.ac.uk/vol1/x.gz", "https://ftp.sra.ebi.ac.uk/vol1/x.gz"),
        ("", ""),
    ],
)
def test_as_url(raw, expected):
    assert _as_url(raw) == expected


def test_parse_filereport_splits_paired_files_with_checksums():
    runs = parse_filereport(FILEREPORT)
    assert len(runs) == 1
    run = runs[0]
    assert run.accession == "SRR1039514"
    assert run.library_strategy == "RNA-Seq"
    assert run.read_count == 40021727

    assert len(run.files) == 2
    first, second = run.files
    assert first.filename == "SRR1039514_1.fastq.gz"
    assert first.url.startswith("https://")
    assert first.size == 2377036173
    # ENA publishes a per-file md5, so checksum verification is real here.
    assert first.md5 == "ae35fe0ce13badacc48c65717e811528"
    assert second.md5 == "9ef4fe59d6378c513f933e24d12f6047"
    assert all(f.kind == "raw" for f in run.files)


def test_parse_filereport_skips_rows_without_a_run():
    assert parse_filereport("run_accession\tfastq_ftp\n\t\n") == []


def test_attach_reads_maps_runs_to_geo_samples(monkeypatch):
    dataset = Dataset(
        metadata=DatasetMetadata(
            accession="GSE52778",
            archive="geo",
            relations={"BioProject": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA229998"},
        ),
        samples=[
            SampleMetadata(
                accession="GSM1275862",
                relations={"BioSample": "https://www.ncbi.nlm.nih.gov/biosample/SAMN02422669"},
            )
        ],
    )
    monkeypatch.setattr(
        "omixfetch.sources.ena.fetch_runs", lambda session, acc: parse_filereport(FILEREPORT)
    )

    runs = attach_reads(dataset, session=None)
    assert len(runs) == 1
    # The run is matched back to its GEO sample through the BioSample id.
    assert runs[0].sample == "GSM1275862"
    assert all(f.accession == "GSM1275862" for f in runs[0].files)
    # Files are added to the dataset manifest.
    assert len(dataset.files) == 2


def test_attach_reads_without_a_project_accession_is_an_error():
    dataset = Dataset(metadata=DatasetMetadata(accession="GSE1", archive="geo"))
    with pytest.raises(ResolveError, match="no BioProject"):
        attach_reads(dataset, session=None)


def test_attach_reads_tolerates_a_project_with_no_public_runs(monkeypatch):
    """Controlled-access series have a project but no ENA runs."""
    dataset = Dataset(
        metadata=DatasetMetadata(
            accession="GSE164378",
            archive="geo",
            relations={"BioProject": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA690251"},
        )
    )
    monkeypatch.setattr("omixfetch.sources.ena.fetch_runs", lambda session, acc: [])
    assert attach_reads(dataset, session=None) == []
    assert dataset.files == []
