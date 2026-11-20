"""Resolve ArrayExpress (E-*) accessions via the EBI BioStudies API.

ArrayExpress now lives inside BioStudies. Study metadata comes from the
submission JSON, per-sample annotation from the MAGE-TAB SDRF, and
study-level annotation from the IDF. File records carry a declared size but
no checksum -- BioStudies does not currently publish md5 sums.
"""

from __future__ import annotations

import csv
import io
from typing import Any

import requests

from ..errors import ResolveError
from ..http import get_json, get_text
from ..models import Dataset, DatasetMetadata, RemoteFile, SampleMetadata

API = "https://www.ebi.ac.uk/biostudies/api/v1/studies"
FILES = "https://www.ebi.ac.uk/biostudies/files"


def file_url(accession: str, path: str) -> str:
    """Return the HTTPS download URL for *path* within *accession*."""
    return f"{FILES}/{accession}/{path.lstrip('/')}"


def _attributes(node: dict[str, Any]) -> dict[str, str]:
    """Flatten a BioStudies ``attributes`` list into a dictionary.

    Repeated names are joined with ``"; "`` so nothing is silently dropped.
    """
    out: dict[str, str] = {}
    for attr in node.get("attributes") or []:
        name, value = attr.get("name"), attr.get("value")
        if not name or value is None:
            continue
        out[name] = f"{out[name]}; {value}" if name in out else str(value)
    return out


def _walk_files(node: Any, collected: list[dict[str, Any]]) -> None:
    """Recursively collect file records; BioStudies nests them in subsections."""
    if isinstance(node, list):
        for item in node:
            _walk_files(item, collected)
    elif isinstance(node, dict):
        entries = node.get("files")
        if entries:
            _flatten_files(entries, collected)
        _walk_files(node.get("subsections") or [], collected)


def _flatten_files(entries: Any, collected: list[dict[str, Any]]) -> None:
    if isinstance(entries, list):
        for item in entries:
            _flatten_files(item, collected)
    elif isinstance(entries, dict):
        collected.append(entries)


def parse_idf(text: str) -> dict[str, list[str]]:
    """Parse a MAGE-TAB IDF into ``{key: [values]}``."""
    out: dict[str, list[str]] = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        cells = line.rstrip("\n").split("\t")
        key = cells[0].strip()
        if not key:
            continue
        values = [c.strip() for c in cells[1:] if c.strip()]
        if values:
            out.setdefault(key, []).extend(values)
    return out


def parse_sdrf(text: str) -> list[dict[str, str]]:
    """Parse a MAGE-TAB SDRF into a list of header-keyed row dictionaries.

    Repeated column headers are suffixed ``#2``, ``#3`` ... so that, for
    example, the several ``Comment[...]`` columns stay distinguishable.
    """
    reader = csv.reader(io.StringIO(text), delimiter="\t")
    rows = [r for r in reader if any(c.strip() for c in r)]
    if not rows:
        return []

    header: list[str] = []
    counts: dict[str, int] = {}
    for col in rows[0]:
        name = col.strip()
        counts[name] = counts.get(name, 0) + 1
        header.append(name if counts[name] == 1 else f"{name}#{counts[name]}")

    return [
        {header[i]: (row[i].strip() if i < len(row) else "") for i in range(len(header))}
        for row in rows[1:]
    ]


def _sdrf_to_samples(rows: list[dict[str, str]]) -> list[SampleMetadata]:
    """Convert SDRF rows into harmonised sample records."""
    samples: list[SampleMetadata] = []
    for row in rows:
        characteristics = {
            key[key.index("[") + 1 : key.rindex("]")]: value
            for key, value in row.items()
            if key.lower().startswith("characteristics") and "[" in key and value
        }
        relations = {
            key[key.index("[") + 1 : key.rindex("]")]: value
            for key, value in row.items()
            if key.lower().startswith("comment") and "[" in key and value
        }
        samples.append(
            SampleMetadata(
                accession=row.get("Source Name", "") or row.get("Assay Name", ""),
                title=row.get("Assay Name") or row.get("Source Name"),
                organism=characteristics.get("organism"),
                source_name=row.get("Source Name"),
                characteristics=characteristics,
                library_strategy=row.get("Comment[LIBRARY_STRATEGY]"),
                library_source=row.get("Comment[LIBRARY_SOURCE]"),
                library_selection=row.get("Comment[LIBRARY_SELECTION]"),
                instrument_model=row.get("Comment[INSTRUMENT_MODEL]"),
                relations=relations,
                raw=dict(row),
            )
        )
    return samples


def resolve_arrayexpress(
    accession: str,
    session: requests.Session,
    *,
    include_samples: bool = True,
    **_: Any,
) -> Dataset:
    """Resolve an ArrayExpress accession into a :class:`Dataset`.

    Args:
        accession: an ArrayExpress accession such as ``E-MTAB-5061``.
        session: HTTP session from :func:`omixfetch.http.make_session`.
        include_samples: download and parse the SDRF for per-sample metadata.
    """
    payload = get_json(session, f"{API}/{accession}")
    if not isinstance(payload, dict) or "section" not in payload:
        raise ResolveError(f"BioStudies returned no usable record for {accession}.")

    top = _attributes(payload)
    section = payload.get("section", {})
    study = _attributes(section)
    study.update({k: v for k, v in top.items() if k not in study})

    records: list[dict[str, Any]] = []
    _walk_files(section, records)

    files: list[RemoteFile] = []
    for record in records:
        path = record.get("path") or record.get("name")
        if not path:
            continue
        attrs = _attributes(record)
        files.append(
            RemoteFile(
                url=file_url(accession, path),
                filename=str(path).rsplit("/", 1)[-1],
                size=record.get("size"),
                kind="metadata"
                if str(path).lower().endswith((".idf.txt", ".sdrf.txt"))
                else "supplementary",
                accession=accession,
                description=attrs.get("Description") or attrs.get("Type"),
            )
        )

    raw_documents: dict[str, str] = {}
    samples: list[SampleMetadata] = []
    protocols: dict[str, str] = {}

    idf_file = next((f for f in files if f.filename.lower().endswith(".idf.txt")), None)
    if idf_file is not None:
        try:
            idf_text = get_text(session, idf_file.url, expect_text=False)
        except ResolveError:
            idf_text = ""
        if idf_text:
            raw_documents[idf_file.filename] = idf_text
            idf = parse_idf(idf_text)
            study.setdefault("Investigation Title", "; ".join(idf.get("Investigation Title", [])))
            names = idf.get("Protocol Name", [])
            descriptions = idf.get("Protocol Description", [])
            protocols = dict(zip(names, descriptions))

    sdrf_file = next(
        (f for f in files if f.filename.lower().endswith(".sdrf.txt")), None
    )
    if include_samples and sdrf_file is not None:
        try:
            sdrf_text = get_text(session, sdrf_file.url, expect_text=False)
        except ResolveError:
            sdrf_text = ""
        if sdrf_text:
            raw_documents[sdrf_file.filename] = sdrf_text
            samples = _sdrf_to_samples(parse_sdrf(sdrf_text))

    organisms = sorted({s.organism for s in samples if s.organism}) or (
        [study["Organism"]] if study.get("Organism") else []
    )

    metadata = DatasetMetadata(
        accession=accession,
        archive="arrayexpress",
        title=study.get("Title") or study.get("Investigation Title"),
        summary=study.get("Description"),
        design=study.get("Study Design") or study.get("Experimental Design"),
        study_type=[study["Study Type"]] if study.get("Study Type") else [],
        organisms=organisms,
        contributors=[a.get("value", "") for a in (section.get("authors") or [])],
        release_date=study.get("ReleaseDate") or study.get("Public Release Date"),
        pubmed_ids=[study[k] for k in ("PubMed ID", "PMID") if study.get(k)],
        protocols=protocols,
        n_samples=len(samples),
        raw={"attributes": study},
    )

    return Dataset(
        metadata=metadata, samples=samples, files=files, raw_documents=raw_documents
    )
