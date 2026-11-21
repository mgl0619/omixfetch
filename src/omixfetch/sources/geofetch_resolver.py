"""Resolve GEO accessions through geofetch instead of the built-in parser.

geofetch (Khoroshevskyi et al., Bioinformatics 2023) maintains its own SOFT
parsing and emits PEP-format metadata. Using it here is optional: install
with ``pip install "omixfetch[geofetch]"``. The built-in resolver stays the
default so the base install remains three dependencies rather than the whole
pepkit stack.

Only metadata resolution is delegated. Transfers still go through omixfetch's
aria2c backend, and raw reads still come from ENA -- geofetch's own raw-data
path shells out to the SRA Toolkit, which pip does not install.
"""

from __future__ import annotations

from typing import Any

import requests

from ..errors import ResolveError
from ..models import Dataset, DatasetMetadata, RemoteFile, SampleMetadata
from .filelist import fetch_filelist

#: Columns geofetch uses for the download URL and size.
URL_COLUMN = "file_url"
SIZE_COLUMN = "file_size"
ACCESSION_COLUMN = "sample_geo_accession"


def geofetch_available() -> bool:
    """Return True if geofetch can be imported."""
    try:
        import geofetch  # noqa: F401
    except ImportError:
        return False
    return True


def _require_geofetch():
    try:
        from geofetch import Geofetcher
    except ImportError as exc:
        raise ResolveError(
            "The geofetch resolver was requested but geofetch is not installed. "
            'Install it with: pip install "omixfetch[geofetch]"'
        ) from exc
    return Geofetcher


def _https(url: str) -> str:
    """geofetch reports GEO URLs with an ftp:// scheme; aria2 wants HTTPS."""
    url = (url or "").strip()
    if url.startswith("ftp://ftp.ncbi.nlm.nih.gov"):
        return "https://ftp.ncbi.nlm.nih.gov" + url[len("ftp://ftp.ncbi.nlm.nih.gov") :]
    return url


def _cell(row: Any, key: str) -> str:
    """Read a cell as a string, mapping pandas NA to an empty string."""
    value = row.get(key) if hasattr(row, "get") else None
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "<na>"} else text


def _size(row: Any) -> int | None:
    raw = _cell(row, SIZE_COLUMN)
    try:
        return int(float(raw)) if raw else None
    except ValueError:
        return None


def _rows(project: Any) -> list[dict[str, Any]]:
    """Return a peppy project's sample table as plain dictionaries."""
    table = getattr(project, "sample_table", None)
    if table is None:
        return []
    return list(table.to_dict(orient="records"))


def _series_metadata(accession: str, project: Any) -> DatasetMetadata:
    """Build study metadata from the series project's config block."""
    config = dict(getattr(project, "config", {}) or {})
    meta = dict(config.get("experiment_metadata") or {})

    def field(name: str) -> str | None:
        return meta.get(f"series_{name}") or None

    def listed(name: str) -> list[str]:
        value = meta.get(f"series_{name}")
        if not value:
            return []
        return [value] if isinstance(value, str) else list(value)

    return DatasetMetadata(
        accession=accession,
        archive="geo",
        title=field("title"),
        summary=field("summary"),
        design=field("overall_design"),
        study_type=listed("type"),
        organisms=listed("sample_organism") or listed("platform_organism"),
        platforms=listed("platform_id"),
        contributors=listed("contributor"),
        submission_date=field("submission_date"),
        release_date=(field("status") or "").replace("Public on ", "") or None,
        last_update_date=field("last_update_date"),
        pubmed_ids=listed("pubmed_id"),
        relations={"raw": field("relation")} if field("relation") else {},
        contact={
            key.replace("series_contact_", ""): str(value)
            for key, value in meta.items()
            if key.startswith("series_contact_")
        },
        raw={"experiment_metadata": meta},
    )


def resolve_geo_via_geofetch(
    accession: str,
    session: requests.Session,
    *,
    include_samples: bool = True,
    include_reads: bool = False,
    prefer: str = "fastq",
    **_: Any,
) -> Dataset:
    """Resolve *accession* using geofetch, returning omixfetch's Dataset.

    Args:
        accession: a GSE accession. geofetch is series-oriented; GSM and GPL
            accessions are not supported by this backend.
        session: HTTP session, used for the ENA step when reads are requested.
        include_samples: include per-sample files and annotation.
        include_reads: additionally resolve raw reads through ENA.

    Raises:
        ResolveError: if geofetch is missing or returns no usable project.
    """
    Geofetcher = _require_geofetch()

    if not accession.upper().startswith("GSE"):
        raise ResolveError(
            f"The geofetch resolver handles GSE series only; {accession} is not one. "
            f"Use --resolver native for GSM and GPL accessions."
        )

    fetcher = Geofetcher(
        processed=True,
        just_metadata=True,
        data_source="all",
        disable_progressbar=True,
    )
    try:
        projects = fetcher.get_projects(
            accession, just_metadata=True, discard_soft=True
        )
    except Exception as exc:  # geofetch raises a variety of types
        raise ResolveError(f"geofetch failed to resolve {accession}: {exc}") from exc

    if not projects:
        raise ResolveError(f"geofetch returned no project for {accession}.")

    series_project = projects.get(f"{accession}_series")
    samples_project = projects.get(f"{accession}_samples")

    metadata = (
        _series_metadata(accession, series_project)
        if series_project is not None
        else DatasetMetadata(accession=accession, archive="geo")
    )

    files: list[RemoteFile] = []
    for row in _rows(series_project) if series_project is not None else []:
        url = _https(_cell(row, URL_COLUMN))
        if url:
            files.append(
                RemoteFile(
                    url=url,
                    filename=url.rsplit("/", 1)[-1],
                    size=_size(row),
                    kind="supplementary",
                    accession=accession,
                )
            )

    samples: list[SampleMetadata] = []
    if include_samples and samples_project is not None:
        # geofetch emits one row per FILE, so rows are grouped back by sample.
        grouped: dict[str, list[dict[str, Any]]] = {}
        for row in _rows(samples_project):
            grouped.setdefault(_cell(row, ACCESSION_COLUMN) or accession, []).append(row)

        for sample_accession, sample_rows in grouped.items():
            first = sample_rows[0]
            sample_files = []
            for row in sample_rows:
                url = _https(_cell(row, URL_COLUMN))
                if url:
                    sample_files.append(
                        RemoteFile(
                            url=url,
                            filename=url.rsplit("/", 1)[-1],
                            size=_size(row),
                            kind="supplementary",
                            accession=sample_accession,
                        )
                    )
            samples.append(
                SampleMetadata(
                    accession=sample_accession,
                    title=_cell(first, "sample_title") or None,
                    organism=_cell(first, "sample_organism_ch1") or None,
                    source_name=_cell(first, "sample_source_name_ch1") or None,
                    characteristics={
                        key.replace("sample_characteristics_", ""): str(value)
                        for key, value in first.items()
                        if key.startswith("sample_characteristics_") and value
                    },
                    library_strategy=_cell(first, "sample_library_strategy") or None,
                    library_source=_cell(first, "sample_library_source") or None,
                    library_selection=_cell(first, "sample_library_selection") or None,
                    instrument_model=_cell(first, "sample_instrument_model") or None,
                    files=sample_files,
                    raw=dict(first),
                )
            )
            files.extend(sample_files)

    metadata.n_samples = len(samples)

    # geofetch does not report a size for series-level files; GEO's own
    # filelist.txt does, and it is one request.
    if any(f.size is None for f in files):
        from .geo import suppl_base_url

        entries = fetch_filelist(session, suppl_base_url(accession, "series"))
        sizes = {e.name: e.size for e in entries if e.size is not None}
        for f in files:
            if f.size is None and f.filename in sizes:
                f.size = sizes[f.filename]

        # filelist.txt covers the RAW archive's members, not necessarily the
        # series-level files, so anything still unknown falls back to a HEAD.
        missing = [f for f in files if f.size is None]
        if missing:
            from .geo import _fill_sizes

            _fill_sizes(session, missing, 8)

    dataset = Dataset(metadata=metadata, samples=samples, files=files)

    if include_reads:
        from .ena import attach_reads

        attach_reads(dataset, session, prefer=prefer)

    return dataset
