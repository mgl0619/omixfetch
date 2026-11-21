"""Write manifests and metadata to disk.

Every download produces, alongside the data, a metadata directory holding the
parsed study and sample annotation plus the unmodified archive records, so a
dataset stays interpretable after the fact.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Sequence
from pathlib import Path

from .models import Dataset, RemoteFile

#: Columns written first in a samples table, when present.
SAMPLE_CORE_COLUMNS = [
    "accession",
    "title",
    "organism",
    "source_name",
    "library_strategy",
    "library_source",
    "library_selection",
    "instrument_model",
]


def write_json(dataset: Dataset, path: Path) -> Path:
    """Write the complete parsed record as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(dataset.to_dict(), indent=2, default=str))
    return path


def write_samples_table(dataset: Dataset, path: Path) -> Path:
    """Write one row per sample, with characteristics flattened into columns.

    The column set is the union across samples, so sparsely annotated samples
    simply carry empty cells rather than being dropped.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    characteristic_keys: list[str] = []
    relation_keys: list[str] = []
    for sample in dataset.samples:
        for key in sample.characteristics:
            if key not in characteristic_keys:
                characteristic_keys.append(key)
        for key in sample.relations:
            if key not in relation_keys:
                relation_keys.append(key)

    header = (
        SAMPLE_CORE_COLUMNS
        + [f"characteristic.{k}" for k in characteristic_keys]
        + [f"relation.{k}" for k in relation_keys]
        + ["n_files"]
    )

    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        for sample in dataset.samples:
            row = [getattr(sample, c, "") or "" for c in SAMPLE_CORE_COLUMNS]
            row += [sample.characteristics.get(k, "") for k in characteristic_keys]
            row += [sample.relations.get(k, "") for k in relation_keys]
            row.append(len(sample.files))
            writer.writerow(row)
    return path


def write_files_table(files: Sequence[RemoteFile], path: Path) -> Path:
    """Write the file manifest as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["accession", "filename", "kind", "size_bytes", "md5", "description", "url"]
        )
        for f in files:
            writer.writerow(
                [
                    f.accession or "",
                    f.filename,
                    f.kind,
                    f.size if f.size is not None else "",
                    f.md5 or "",
                    f.description or "",
                    f.url,
                ]
            )
    return path


def write_runs_table(dataset: Dataset, path: Path) -> Path:
    """Write one row per sequencing run resolved from ENA."""
    path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "run_accession",
        "geo_sample",
        "experiment_accession",
        "sample_accession",
        "study_accession",
        "library_strategy",
        "library_source",
        "library_selection",
        "instrument_model",
        "read_count",
        "base_count",
        "n_files",
        "total_bytes",
        "urls",
        "md5s",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        for run in dataset.runs:
            writer.writerow(
                [
                    run.accession,
                    run.sample or "",
                    run.experiment_accession or "",
                    run.sample_accession or "",
                    run.study_accession or "",
                    run.library_strategy or "",
                    run.library_source or "",
                    run.library_selection or "",
                    run.instrument_model or "",
                    run.read_count if run.read_count is not None else "",
                    run.base_count if run.base_count is not None else "",
                    len(run.files),
                    sum(f.size or 0 for f in run.files),
                    ";".join(f.url for f in run.files),
                    ";".join(f.md5 or "" for f in run.files),
                ]
            )
    return path


def write_raw_documents(dataset: Dataset, directory: Path) -> list[Path]:
    """Write the unmodified archive records (SOFT, IDF, SDRF)."""
    directory.mkdir(parents=True, exist_ok=True)
    written = []
    for name, text in dataset.raw_documents.items():
        path = directory / name
        path.write_text(text)
        written.append(path)
    return written


def write_metadata(dataset: Dataset, directory: Path) -> dict[str, Path]:
    """Write every metadata artefact for *dataset* into *directory*.

    Returns:
        Mapping of artefact name to the path written.
    """
    directory = Path(directory)
    acc = dataset.accession
    written = {
        "json": write_json(dataset, directory / f"{acc}.metadata.json"),
        "files": write_files_table(dataset.files, directory / f"{acc}.files.tsv"),
    }
    if dataset.samples:
        written["samples"] = write_samples_table(
            dataset, directory / f"{acc}.samples.tsv"
        )
    if dataset.runs:
        written["runs"] = write_runs_table(dataset, directory / f"{acc}.runs.tsv")
    write_raw_documents(dataset, directory)
    return written


def total_size(files: Sequence[RemoteFile]) -> int | None:
    """Sum known file sizes, or return None when no size is known.

    Returning None rather than 0 keeps "sizes were never queried" distinct
    from "these files are empty".
    """
    known = [f.size for f in files if f.size is not None]
    return sum(known) if known else None


def human_bytes(n: int | None) -> str:
    """Format a byte count for display."""
    if n is None:
        return "unknown"
    size = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.1f} {unit}" if unit != "B" else f"{int(size)} B"
        size /= 1024
    return f"{size:.1f} TB"
