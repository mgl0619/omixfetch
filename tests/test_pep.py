"""Native PEP export."""

import csv

import yaml

from omixfetch.models import Dataset, DatasetMetadata, RemoteFile, SampleMetadata
from omixfetch.pep import PEP_VERSION, write_pep


def _dataset():
    return Dataset(
        metadata=DatasetMetadata(
            accession="GSE000001",
            archive="geo",
            title="A test series",
            organisms=["Homo sapiens"],
            pubmed_ids=["33456401"],
            n_samples=2,
        ),
        samples=[
            SampleMetadata(
                accession="GSM000001",
                title="Control",
                organism="Homo sapiens",
                characteristics={"cell type": "PBMC"},
                relations={"BioSample": "SAMN1"},
                files=[RemoteFile(url="https://x/a.gz", filename="a.gz")],
            ),
            SampleMetadata(
                accession="GSM000002",
                title="Treated",
                characteristics={"cell type": "PBMC", "timepoint": "24h"},
            ),
        ],
    )


def test_write_pep_emits_config_and_table(tmp_path):
    written = write_pep(_dataset(), tmp_path)
    assert set(written) == {"config", "sample_table"}
    assert written["config"].exists() and written["sample_table"].exists()


def test_config_declares_version_and_points_at_the_table(tmp_path):
    written = write_pep(_dataset(), tmp_path)
    config = yaml.safe_load(written["config"].read_text())

    assert config["pep_version"] == PEP_VERSION
    # peppy reads "name"; without it the project is named after its directory.
    assert config["name"] == "GSE000001"
    assert config["project_name"] == "GSE000001"
    assert config["sample_table"] == "GSE000001_samples.csv"
    assert (written["config"].parent / config["sample_table"]).exists()
    assert config["experiment_metadata"]["title"] == "A test series"


def test_config_omits_empty_metadata(tmp_path):
    written = write_pep(_dataset(), tmp_path)
    meta = yaml.safe_load(written["config"].read_text())["experiment_metadata"]
    # Fields with no value are left out rather than written as nulls.
    assert "summary" not in meta
    assert "design" not in meta


def test_sample_table_has_sample_name_and_flattened_characteristics(tmp_path):
    written = write_pep(_dataset(), tmp_path)
    rows = list(csv.DictReader(written["sample_table"].open()))

    assert len(rows) == 2
    # sample_name is required by the PEP spec.
    assert rows[0]["sample_name"] == "GSM000001"
    # Column names are normalised into identifiers.
    assert rows[0]["char_cell_type"] == "PBMC"
    assert rows[0]["char_timepoint"] == ""
    assert rows[1]["char_timepoint"] == "24h"
    assert rows[0]["rel_biosample"] == "SAMN1"
    assert rows[0]["n_files"] == "1"
    assert rows[0]["file_urls"] == "https://x/a.gz"


def test_unicode_is_written_readably(tmp_path):
    dataset = _dataset()
    dataset.metadata.summary = "\u2018weighted-nearest neighbor\u2019 analysis"
    written = write_pep(dataset, tmp_path)
    # allow_unicode keeps the text legible instead of \uXXXX escapes.
    assert "weighted-nearest neighbor" in written["config"].read_text()
