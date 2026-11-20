"""Resolve NCBI GEO accessions into metadata and downloadable files.

File enumeration goes through the SOFT record rather than an FTP directory
listing: ``https://ftp.ncbi.nlm.nih.gov/geo/.../suppl/`` returns HTTP 403 for
directory requests, while the per-record ``supplementary_file`` fields are
authoritative and come back in a single query.
"""

from __future__ import annotations

import re
from concurrent.futures import ThreadPoolExecutor

import requests

from ..accession import geo_ftp_stub, geo_kind
from ..errors import ResolveError
from ..http import get_text, head_size
from ..models import Dataset, DatasetMetadata, RemoteFile, SampleMetadata

ACC_CGI = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
FTP_HTTPS = "https://ftp.ncbi.nlm.nih.gov"

#: ``!Series_title = Some title`` -> ("Series_title", "Some title")
FIELD_RE = re.compile(r"^!([A-Za-z0-9_]+)\s*=\s*(.*)$")
ENTITY_RE = re.compile(r"^\^([A-Z]+)\s*=\s*(.*)$")


def parse_soft(text: str) -> list[tuple[str, str, dict[str, list[str]]]]:
    """Parse SOFT text into ``(entity_type, accession, fields)`` triples.

    Repeated keys are preserved as lists, because GEO uses repetition for
    multi-valued fields such as ``Sample_characteristics_ch1``.
    """
    records: list[tuple[str, str, dict[str, list[str]]]] = []
    entity_type = accession = None
    fields: dict[str, list[str]] = {}

    for line in text.splitlines():
        entity = ENTITY_RE.match(line)
        if entity:
            if entity_type is not None:
                records.append((entity_type, accession or "", fields))
            entity_type, accession, fields = entity.group(1), entity.group(2).strip(), {}
            continue
        field = FIELD_RE.match(line)
        if field and entity_type is not None:
            fields.setdefault(field.group(1), []).append(field.group(2).strip())

    if entity_type is not None:
        records.append((entity_type, accession or "", fields))
    return records


def _first(fields: dict[str, list[str]], key: str) -> str | None:
    values = fields.get(key)
    return values[0] if values else None


def _split_pairs(values: list[str]) -> dict[str, str]:
    """Split ``"cell type: PBMC"`` style entries into a dictionary."""
    out: dict[str, str] = {}
    for value in values:
        key, sep, rest = value.partition(":")
        if sep:
            out[key.strip()] = rest.strip()
        else:
            out.setdefault("annotation", value.strip())
    return out


def _ftp_to_https(url: str) -> str:
    """Rewrite GEO's ``ftp://`` URLs to HTTPS (the sandboxed FTP port is closed
    in many environments, and HTTPS supports range requests)."""
    return re.sub(r"^ftp://ftp\.ncbi\.nlm\.nih\.gov", FTP_HTTPS, url.strip())


def _supplementary_files(
    fields: dict[str, list[str]], prefix: str, accession: str
) -> list[RemoteFile]:
    """Collect every ``*_supplementary_file*`` entry from a SOFT record."""
    files: list[RemoteFile] = []
    for key, values in fields.items():
        if not key.startswith(f"{prefix}_supplementary_file"):
            continue
        for value in values:
            if not value or value.upper() == "NONE":
                continue
            url = _ftp_to_https(value)
            if not url.lower().startswith("http"):
                continue
            files.append(
                RemoteFile(
                    url=url,
                    filename=url.rsplit("/", 1)[-1],
                    kind="supplementary",
                    accession=accession,
                )
            )
    return files


def series_matrix_url(accession: str) -> str:
    """Return the canonical series-matrix URL for a GSE accession."""
    stub = geo_ftp_stub(accession)
    return (
        f"{FTP_HTTPS}/geo/series/{stub}/{accession}/matrix/"
        f"{accession}_series_matrix.txt.gz"
    )


def _sample_metadata(accession: str, fields: dict[str, list[str]]) -> SampleMetadata:
    relations = _split_pairs(fields.get("Sample_relation", []))
    return SampleMetadata(
        accession=accession,
        title=_first(fields, "Sample_title"),
        organism=_first(fields, "Sample_organism_ch1"),
        source_name=_first(fields, "Sample_source_name_ch1"),
        characteristics=_split_pairs(fields.get("Sample_characteristics_ch1", [])),
        library_strategy=_first(fields, "Sample_library_strategy"),
        library_source=_first(fields, "Sample_library_source"),
        library_selection=_first(fields, "Sample_library_selection"),
        instrument_model=_first(fields, "Sample_instrument_model"),
        relations=relations,
        files=_supplementary_files(fields, "Sample", accession),
        raw=dict(fields),
    )


def _series_metadata(accession: str, fields: dict[str, list[str]]) -> DatasetMetadata:
    contact = {
        key.replace("Series_contact_", ""): values[0]
        for key, values in fields.items()
        if key.startswith("Series_contact_") and values
    }
    return DatasetMetadata(
        accession=accession,
        archive="geo",
        title=_first(fields, "Series_title"),
        summary=" ".join(fields.get("Series_summary", [])) or None,
        design=_first(fields, "Series_overall_design"),
        study_type=fields.get("Series_type", []),
        organisms=sorted(set(fields.get("Series_sample_organism", []))),
        platforms=fields.get("Series_platform_id", []),
        contributors=fields.get("Series_contributor", []),
        submission_date=_first(fields, "Series_submission_date"),
        release_date=(_first(fields, "Series_status") or "").replace("Public on ", "")
        or None,
        last_update_date=_first(fields, "Series_last_update_date"),
        pubmed_ids=fields.get("Series_pubmed_id", []),
        relations=_split_pairs(fields.get("Series_relation", [])),
        contact=contact,
        n_samples=len(fields.get("Series_sample_id", [])),
        raw=dict(fields),
    )


def resolve_geo(
    accession: str,
    session: requests.Session,
    *,
    include_samples: bool = True,
    include_matrix: bool = False,
    fetch_sizes: bool = True,
    size_workers: int = 8,
) -> Dataset:
    """Resolve a GEO accession into a :class:`~omixfetch.models.Dataset`.

    Args:
        accession: a GSE, GSM or GPL accession.
        session: HTTP session from :func:`omixfetch.http.make_session`.
        include_samples: for a series, also fetch every sample's metadata and
            per-sample supplementary files. This is one extra request total,
            not one per sample.
        include_matrix: add the series-matrix file to the download list.
        fetch_sizes: issue HEAD requests to learn file sizes, which GEO does
            not publish in SOFT.
        size_workers: concurrency for those HEAD requests.
    """
    kind = geo_kind(accession)

    self_soft = get_text(
        session,
        ACC_CGI,
        params={"acc": accession, "targ": "self", "form": "text", "view": "brief"},
    )
    records = parse_soft(self_soft)
    if not records:
        raise ResolveError(f"GEO returned no SOFT record for {accession}.")

    _, _, fields = records[0]
    raw_documents = {f"{accession}.soft.txt": self_soft}
    samples: list[SampleMetadata] = []
    files: list[RemoteFile] = []

    if kind == "series":
        metadata = _series_metadata(accession, fields)
        files.extend(_supplementary_files(fields, "Series", accession))

        if include_samples and metadata.n_samples:
            sample_soft = get_text(
                session,
                ACC_CGI,
                params={
                    "acc": accession,
                    "targ": "gsm",
                    "form": "text",
                    "view": "brief",
                },
            )
            raw_documents[f"{accession}.samples.soft.txt"] = sample_soft
            for entity, sample_acc, sample_fields in parse_soft(sample_soft):
                if entity == "SAMPLE":
                    sample = _sample_metadata(sample_acc, sample_fields)
                    samples.append(sample)
                    files.extend(sample.files)

        if include_matrix:
            files.append(
                RemoteFile(
                    url=series_matrix_url(accession),
                    filename=f"{accession}_series_matrix.txt.gz",
                    kind="matrix",
                    accession=accession,
                )
            )
    else:
        prefix = "Sample" if kind == "samples" else "Platform"
        title = _first(fields, f"{prefix}_title")
        metadata = DatasetMetadata(
            accession=accession,
            archive="geo",
            title=title,
            organisms=fields.get(f"{prefix}_organism_ch1", [])
            or fields.get("Platform_organism", []),
            last_update_date=_first(fields, f"{prefix}_last_update_date"),
            raw=dict(fields),
        )
        files.extend(_supplementary_files(fields, prefix, accession))
        if kind == "samples":
            samples.append(_sample_metadata(accession, fields))
            metadata.n_samples = 1

    files = _dedupe(files)
    if fetch_sizes and files:
        _fill_sizes(session, files, size_workers)

    return Dataset(
        metadata=metadata, samples=samples, files=files, raw_documents=raw_documents
    )


def _dedupe(files: list[RemoteFile]) -> list[RemoteFile]:
    """Drop duplicate URLs, preserving order."""
    seen: set[str] = set()
    out: list[RemoteFile] = []
    for f in files:
        if f.url not in seen:
            seen.add(f.url)
            out.append(f)
    return out


def _fill_sizes(
    session: requests.Session, files: list[RemoteFile], workers: int
) -> None:
    """Populate ``size`` on each file via concurrent HEAD requests."""
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        for f, size in zip(files, pool.map(lambda x: head_size(session, x.url), files)):
            f.size = size
