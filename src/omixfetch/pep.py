"""Export a resolved dataset as a PEP (Portable Encapsulated Project).

PEP is a two-file sample-metadata standard from pepkit: a YAML config naming
the project and pointing at a sample table, plus the table itself. Tools in
that ecosystem -- peppy, looper, eido, pepr -- consume the pair directly.

Writing it here rather than through ``peppy`` is deliberate: peppy pulls in
numpy, pandas and the rest of the pepkit stack, while the format itself is a
small YAML document and a CSV.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import yaml

from .models import Dataset

#: PEP specification version the exported config declares.
PEP_VERSION = "2.1.0"


def _project_metadata(dataset: Dataset) -> dict[str, Any]:
    """Collect study-level annotation for the config's metadata block."""
    m = dataset.metadata
    fields = {
        "accession": m.accession,
        "archive": m.archive,
        "title": m.title,
        "summary": m.summary,
        "design": m.design,
        "study_type": m.study_type,
        "organisms": m.organisms,
        "platforms": m.platforms,
        "contributors": m.contributors,
        "submission_date": m.submission_date,
        "release_date": m.release_date,
        "last_update_date": m.last_update_date,
        "pubmed_ids": m.pubmed_ids,
        "relations": m.relations,
        "n_samples": m.n_samples,
    }
    return {k: v for k, v in fields.items() if v not in (None, "", [], {})}


def _sample_rows(dataset: Dataset) -> tuple[list[str], list[dict[str, Any]]]:
    """Build the sample table: one row per sample, characteristics as columns.

    ``sample_name`` is required by the PEP spec and is taken from the archive
    accession, which is already unique within a dataset.
    """
    characteristic_keys: list[str] = []
    relation_keys: list[str] = []
    for sample in dataset.samples:
        for key in sample.characteristics:
            if key not in characteristic_keys:
                characteristic_keys.append(key)
        for key in sample.relations:
            if key not in relation_keys:
                relation_keys.append(key)

    def column(prefix: str, key: str) -> str:
        # PEP consumers treat column names as identifiers; spaces and dots
        # make them awkward to reference in schemas and pipelines.
        safe = key.strip().lower().replace(" ", "_").replace(".", "_")
        return f"{prefix}_{safe}"

    header = [
        "sample_name",
        "accession",
        "title",
        "organism",
        "source_name",
        "library_strategy",
        "library_source",
        "library_selection",
        "instrument_model",
    ]
    header += [column("char", k) for k in characteristic_keys]
    header += [column("rel", k) for k in relation_keys]
    header += ["n_files", "file_urls"]

    rows: list[dict[str, Any]] = []
    for sample in dataset.samples:
        row: dict[str, Any] = {
            "sample_name": sample.accession,
            "accession": sample.accession,
            "title": sample.title or "",
            "organism": sample.organism or "",
            "source_name": sample.source_name or "",
            "library_strategy": sample.library_strategy or "",
            "library_source": sample.library_source or "",
            "library_selection": sample.library_selection or "",
            "instrument_model": sample.instrument_model or "",
        }
        for key in characteristic_keys:
            row[column("char", key)] = sample.characteristics.get(key, "")
        for key in relation_keys:
            row[column("rel", key)] = sample.relations.get(key, "")
        row["n_files"] = len(sample.files)
        row["file_urls"] = ";".join(f.url for f in sample.files)
        rows.append(row)

    return header, rows


def write_pep(dataset: Dataset, directory: Path) -> dict[str, Path]:
    """Write ``<ACC>_config.yaml`` and ``<ACC>_samples.csv`` into *directory*.

    Returns:
        Mapping of ``"config"`` and ``"sample_table"`` to the paths written.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    accession = dataset.accession

    table_name = f"{accession}_samples.csv"
    header, rows = _sample_rows(dataset)
    table_path = directory / table_name
    with open(table_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)

    config = {
        "pep_version": PEP_VERSION,
        # PEP 2.1.0 reads the project name from "name"; peppy otherwise falls
        # back to the containing directory name. "project_name" is kept too
        # because geofetch-produced PEPs use it and tooling reads both.
        "name": accession,
        "project_name": accession,
        "sample_table": table_name,
        "experiment_metadata": _project_metadata(dataset),
    }
    config_path = directory / f"{accession}_config.yaml"
    with open(config_path, "w") as handle:
        yaml.safe_dump(
            config,
            handle,
            default_flow_style=False,
            sort_keys=False,
            allow_unicode=True,
        )

    return {"config": config_path, "sample_table": table_path}
