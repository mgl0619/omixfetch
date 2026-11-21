"""Dataclasses describing resolved datasets."""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Any


@dataclass
class RemoteFile:
    """One downloadable file belonging to a dataset."""

    url: str
    filename: str
    #: Declared size in bytes, or ``None`` when the archive does not publish one.
    size: int | None = None
    #: ``"supplementary"``, ``"metadata"``, ``"matrix"`` or ``"raw"``.
    kind: str = "supplementary"
    #: Accession this file hangs off (the series, or an individual sample).
    accession: str | None = None
    #: md5 as published by the archive. GEO does not publish checksums, so
    #: this is ``None`` for every GEO file; the field exists for archives and
    #: backends that do supply one.
    md5: str | None = None
    #: Free-text description from the archive, when available.
    description: str | None = None
    #: Path relative to the download directory, set by the layout step.
    #: ``None`` means the file lands at the top level under its own name.
    relpath: str | None = None

    @property
    def target(self) -> str:
        """Relative path this file should be written to."""
        return self.relpath or self.filename

    def to_dict(self) -> dict[str, Any]:
        return dataclasses.asdict(self)


@dataclass
class SampleMetadata:
    """Per-sample annotation, harmonised across archives."""

    accession: str
    title: str | None = None
    organism: str | None = None
    source_name: str | None = None
    #: Free-form key/value annotation parsed from GEO characteristics fields.
    characteristics: dict[str, str] = field(default_factory=dict)
    library_strategy: str | None = None
    library_source: str | None = None
    library_selection: str | None = None
    instrument_model: str | None = None
    #: Linked archive identifiers, e.g. ENA run or BioSample accessions.
    relations: dict[str, str] = field(default_factory=dict)
    files: list[RemoteFile] = field(default_factory=list)
    #: Everything the archive returned, unmodified.
    raw: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return dataclasses.asdict(self)


@dataclass
class DatasetMetadata:
    """Study-level annotation, harmonised across archives."""

    accession: str
    archive: str
    title: str | None = None
    summary: str | None = None
    design: str | None = None
    study_type: list[str] = field(default_factory=list)
    organisms: list[str] = field(default_factory=list)
    platforms: list[str] = field(default_factory=list)
    contributors: list[str] = field(default_factory=list)
    submission_date: str | None = None
    release_date: str | None = None
    last_update_date: str | None = None
    pubmed_ids: list[str] = field(default_factory=list)
    relations: dict[str, str] = field(default_factory=dict)
    contact: dict[str, str] = field(default_factory=dict)
    protocols: dict[str, str] = field(default_factory=dict)
    n_samples: int = 0
    #: Which resolver produced this record ("native" or "geofetch"), so a bug
    #: report says which code path ran.
    resolved_by: str = "native"
    raw: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return dataclasses.asdict(self)


@dataclass
class Dataset:
    """A fully resolved dataset: metadata, samples and downloadable files."""

    metadata: DatasetMetadata
    samples: list[SampleMetadata] = field(default_factory=list)
    files: list[RemoteFile] = field(default_factory=list)
    #: Raw metadata documents (SOFT, IDF, SDRF) keyed by filename, kept so the
    #: exact archive record can be written to disk alongside the parsed form.
    raw_documents: dict[str, str] = field(default_factory=dict)
    #: Sequencing runs resolved from ENA. Typed as Any to avoid a circular
    #: import; entries are omixfetch.sources.ena.RunMetadata.
    runs: list[Any] = field(default_factory=list)

    @property
    def accession(self) -> str:
        return self.metadata.accession

    @property
    def total_bytes(self) -> int | None:
        """Sum of declared file sizes, or None when no size is known."""
        known = [f.size for f in self.files if f.size is not None]
        return sum(known) if known else None

    def to_dict(self) -> dict[str, Any]:
        return {
            "metadata": self.metadata.to_dict(),
            "samples": [s.to_dict() for s in self.samples],
            "files": [f.to_dict() for f in self.files],
            "runs": [r.to_dict() for r in self.runs],
        }
