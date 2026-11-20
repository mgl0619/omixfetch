"""Metadata and manifest writers."""

import csv
import json

from omixfetch.manifest import human_bytes, total_size, write_files_table, write_metadata
from omixfetch.models import Dataset, DatasetMetadata, RemoteFile, SampleMetadata


def _dataset():
    return Dataset(
        metadata=DatasetMetadata(
            accession="GSE000001", archive="geo", title="A test", n_samples=2
        ),
        samples=[
            SampleMetadata(
                accession="GSM000001",
                title="s1",
                characteristics={"cell type": "PBMC"},
                relations={"BioSample": "SAMN1"},
            ),
            SampleMetadata(
                accession="GSM000002",
                title="s2",
                characteristics={"cell type": "PBMC", "timepoint": "24h"},
            ),
        ],
        files=[RemoteFile(url="https://x/a.tar", filename="a.tar", size=100)],
        raw_documents={"GSE000001.soft.txt": "^SERIES = GSE000001\n"},
    )


def test_write_metadata_emits_every_artefact(tmp_path):
    written = write_metadata(_dataset(), tmp_path)
    assert set(written) == {"json", "files", "samples"}
    assert (tmp_path / "GSE000001.metadata.json").exists()
    assert (tmp_path / "GSE000001.samples.tsv").exists()
    assert (tmp_path / "GSE000001.files.tsv").exists()
    # The unmodified archive record is kept alongside the parsed form.
    assert (tmp_path / "GSE000001.soft.txt").read_text().startswith("^SERIES")


def test_samples_table_unions_sparse_characteristics(tmp_path):
    write_metadata(_dataset(), tmp_path)
    with open(tmp_path / "GSE000001.samples.tsv") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 2
    assert "characteristic.timepoint" in rows[0]
    # Sample 1 lacks timepoint, so the cell is empty rather than missing.
    assert rows[0]["characteristic.timepoint"] == ""
    assert rows[1]["characteristic.timepoint"] == "24h"
    assert rows[0]["relation.BioSample"] == "SAMN1"


def test_metadata_json_roundtrips(tmp_path):
    write_metadata(_dataset(), tmp_path)
    payload = json.loads((tmp_path / "GSE000001.metadata.json").read_text())
    assert payload["metadata"]["accession"] == "GSE000001"
    assert len(payload["samples"]) == 2
    assert payload["files"][0]["filename"] == "a.tar"


def test_files_table_handles_unknown_size(tmp_path):
    path = write_files_table(
        [RemoteFile(url="https://x/b.gz", filename="b.gz", size=None)],
        tmp_path / "f.tsv",
    )
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    assert rows[0]["size_bytes"] == ""


def test_human_bytes():
    assert human_bytes(None) == "unknown"
    assert human_bytes(512) == "512 B"
    assert human_bytes(1536) == "1.5 KB"
    assert human_bytes(5 * (1 << 30)) == "5.0 GB"


def test_total_size_distinguishes_unknown_from_empty():
    """A manifest whose sizes were never queried must not report 0 B."""
    assert total_size([]) is None
    assert total_size([RemoteFile(url="u", filename="a", size=None)]) is None
    assert total_size(
        [RemoteFile(url="u", filename="a", size=None), RemoteFile(url="v", filename="b", size=7)]
    ) == 7


def test_dataset_total_bytes_is_none_when_unknown():
    dataset = Dataset(
        metadata=DatasetMetadata(accession="GSE1", archive="geo"),
        files=[RemoteFile(url="u", filename="a", size=None)],
    )
    assert dataset.total_bytes is None
