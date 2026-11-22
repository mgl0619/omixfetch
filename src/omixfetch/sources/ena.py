"""Resolve raw sequencing reads through the EBI ENA portal API.

ENA mirrors SRA and serves FASTQ files as ordinary HTTPS objects with range
support, so raw reads can be fetched with the same aria2c machinery as
processed files -- no SRA Toolkit, and no ``.sra`` to convert. ENA also
publishes a per-file md5, which makes checksum verification real rather than
nominal: GEO publishes no checksums at all.

Not every GEO series has reads here. Human data under controlled access
(dbGaP) has a BioProject but no public runs, and ENA then returns an empty
result. That is a permissions boundary, not an error to retry.
"""

from __future__ import annotations

import csv
import dataclasses
import io
import re
from dataclasses import dataclass, field
from typing import Any

import requests

from ..errors import ResolveError
from ..http import get_text
from ..models import Dataset, RemoteFile, SampleMetadata

PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

#: Fields requested from the portal, in the order they are returned.
RUN_FIELDS = [
    "run_accession",
    "experiment_accession",
    "sample_accession",
    "secondary_sample_accession",
    "study_accession",
    "library_strategy",
    "library_source",
    "library_selection",
    "instrument_model",
    "read_count",
    "base_count",
    "fastq_ftp",
    "fastq_aspera",
    "fastq_bytes",
    "fastq_md5",
    "submitted_ftp",
    "submitted_aspera",
    "submitted_bytes",
    "submitted_md5",
]

#: Study-level accessions ENA accepts for a read_run query.
PROJECT_RE = re.compile(r"\b(PRJ[A-Z]{2}\d+)\b")
STUDY_RE = re.compile(r"\b([SED]RP\d+)\b")
BIOSAMPLE_RE = re.compile(r"\b(SAM[ND][A-Z]?\d+|SAME[A]?\d+)\b")


@dataclass
class RunMetadata:
    """One sequencing run as described by ENA."""

    accession: str
    experiment_accession: str | None = None
    sample_accession: str | None = None
    secondary_sample_accession: str | None = None
    study_accession: str | None = None
    library_strategy: str | None = None
    library_source: str | None = None
    library_selection: str | None = None
    instrument_model: str | None = None
    read_count: int | None = None
    base_count: int | None = None
    #: GEO sample this run belongs to, when it could be matched.
    sample: str | None = None
    files: list[RemoteFile] = field(default_factory=list)
    raw: dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return dataclasses.asdict(self)


def project_accessions(relations: dict[str, str]) -> list[str]:
    """Extract ENA-queryable study accessions from GEO relation fields.

    GEO records these as URLs -- ``BioProject: https://.../PRJNA229998`` and
    ``SRA: https://...?term=SRP033351``. Both forms are accepted by the
    portal; BioProject accessions are listed first because they are present
    more often.
    """
    projects: list[str] = []
    studies: list[str] = []
    for value in relations.values():
        projects += PROJECT_RE.findall(value or "")
        studies += STUDY_RE.findall(value or "")

    seen: set[str] = set()
    ordered: list[str] = []
    for acc in projects + studies:
        if acc not in seen:
            seen.add(acc)
            ordered.append(acc)
    return ordered


def _as_int(value: str | None) -> int | None:
    """Parse an integer column, tolerating blanks the portal leaves empty."""
    value = (value or "").strip()
    return int(value) if value.isdigit() else None


def _as_url(path: str) -> str:
    """ENA returns bare host paths; make them absolute HTTPS URLs."""
    path = path.strip()
    if not path:
        return ""
    if path.startswith(("http://", "https://")):
        return path
    if path.startswith("ftp://"):
        return "https://" + path[len("ftp://") :]
    return "https://" + path


#: Column groups for each file flavour the portal serves.
FASTQ_COLUMNS = ("fastq_ftp", "fastq_bytes", "fastq_md5", "fastq_aspera")
SUBMITTED_COLUMNS = (
    "submitted_ftp",
    "submitted_bytes",
    "submitted_md5",
    "submitted_aspera",
)


def _row_files(
    row: dict[str, str], run_accession: str, prefer: str = "fastq"
) -> list[RemoteFile]:
    """Build file records from a portal row.

    Args:
        prefer: ``"fastq"`` for ENA-generated FASTQ, or ``"submitted"`` for the
            files the submitter uploaded (often BAM or CRAM). Either way the
            other flavour is used when the preferred one is absent, so a run
            with only BAMs still resolves.

    The portal packs multi-file runs into semicolon-separated lists that run
    in parallel across the URL, size and md5 columns.
    """
    order = (
        (SUBMITTED_COLUMNS, FASTQ_COLUMNS)
        if prefer == "submitted"
        else (FASTQ_COLUMNS, SUBMITTED_COLUMNS)
    )
    for url_key, size_key, md5_key, aspera_key in order:
        raw_urls = [u for u in (row.get(url_key) or "").split(";") if u.strip()]
        if not raw_urls:
            continue
        sizes = (row.get(size_key) or "").split(";")
        md5s = (row.get(md5_key) or "").split(";")
        aspera = [a for a in (row.get(aspera_key) or "").split(";") if a.strip()]
        files: list[RemoteFile] = []
        for index, raw_url in enumerate(raw_urls):
            url = _as_url(raw_url)
            size: int | None = None
            if index < len(sizes) and sizes[index].strip().isdigit():
                size = int(sizes[index].strip())
            md5 = md5s[index].strip() if index < len(md5s) else ""
            md5 = md5 or None
            files.append(
                RemoteFile(
                    url=url,
                    filename=url.rsplit("/", 1)[-1],
                    size=size,
                    kind="raw",
                    accession=run_accession,
                    md5=md5,
                    aspera=aspera[index].strip() if index < len(aspera) else None,
                )
            )
        return files
    return []


def parse_filereport(text: str, prefer: str = "fastq") -> list[RunMetadata]:
    """Parse a TSV filereport response into run records.

    Args:
        prefer: ``"fastq"`` or ``"submitted"`` -- see :func:`_row_files`.
    """
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    runs: list[RunMetadata] = []
    for row in rows:
        accession = (row.get("run_accession") or "").strip()
        if not accession:
            continue
        runs.append(
            RunMetadata(
                accession=accession,
                experiment_accession=row.get("experiment_accession") or None,
                sample_accession=row.get("sample_accession") or None,
                secondary_sample_accession=row.get("secondary_sample_accession") or None,
                study_accession=row.get("study_accession") or None,
                library_strategy=row.get("library_strategy") or None,
                library_source=row.get("library_source") or None,
                library_selection=row.get("library_selection") or None,
                instrument_model=row.get("instrument_model") or None,
                read_count=_as_int(row.get("read_count")),
                base_count=_as_int(row.get("base_count")),
                files=_row_files(row, accession, prefer=prefer),
                raw={k: v for k, v in row.items() if v},
            )
        )
    return runs


def fetch_runs(
    session: requests.Session, accession: str, prefer: str = "fastq"
) -> list[RunMetadata]:
    """Query the ENA portal for every public run under *accession*.

    Args:
        accession: a BioProject (``PRJNA...``) or study (``SRP...``) accession.

    Returns:
        Run records, which may be empty when the data is not public in ENA.
    """
    text = get_text(
        session,
        PORTAL,
        expect_text=False,
        params={
            "accession": accession,
            "result": "read_run",
            "fields": ",".join(RUN_FIELDS),
            "format": "tsv",
        },
    )
    return parse_filereport(text, prefer=prefer)


def _biosample_to_sample(samples: list[SampleMetadata]) -> dict[str, str]:
    """Map BioSample accessions mentioned in GEO relations to GEO sample ids."""
    mapping: dict[str, str] = {}
    for sample in samples:
        for value in sample.relations.values():
            for match in BIOSAMPLE_RE.findall(value or ""):
                mapping[match] = sample.accession
    return mapping


def attach_reads(
    dataset: Dataset,
    session: requests.Session,
    prefer: str = "fastq",
    source: str = "ena",
    sra_flavour: str = "sra",
) -> list[RunMetadata]:
    """Resolve raw reads for *dataset* and attach them in place.

    Args:
        prefer: for the ENA source, ``"fastq"`` or ``"submitted"``.
        source: where the bytes come from. ``"ena"`` uses ENA's own file
            locations; ``"ncbi"`` keeps ENA for run discovery and sample
            mapping but takes download locations from NCBI's SDL service,
            which is often considerably faster and can serve the smaller
            ``.sra`` archive instead of a large submitted BAM.
        sra_flavour: with ``source="ncbi"``, ``"sra"`` for the normalised
            archive or ``"original"`` for the submitter's upload.

    Returns:
        The runs found. An empty list means the series links to a project but
        ENA serves no public runs for it -- usually controlled access.

    Raises:
        ResolveError: if the dataset carries no project accession at all.
    """
    accessions = project_accessions(dataset.metadata.relations)
    if not accessions:
        raise ResolveError(
            f"{dataset.accession} records no BioProject or SRA study accession, so "
            f"there is no raw-read project to look up. Series with processed data "
            f"only are the usual reason."
        )

    runs: list[RunMetadata] = []
    for accession in accessions:
        runs = fetch_runs(session, accession, prefer=prefer)
        if runs:
            break

    by_biosample = _biosample_to_sample(dataset.samples)
    for run in runs:
        for key in (run.sample_accession, run.secondary_sample_accession):
            if key and key in by_biosample:
                run.sample = by_biosample[key]
                break

    if source == "ncbi":
        _use_ncbi_locations(runs, session, sra_flavour)

    for run in runs:
        for f in run.files:
            f.accession = run.sample or run.accession

    dataset.runs = runs
    dataset.files.extend(f for run in runs for f in run.files)
    return runs


def _use_ncbi_locations(
    runs: list[RunMetadata], session: requests.Session, flavour: str
) -> None:
    """Replace each run's files with NCBI SDL locations, in place.

    Runs SDL has nothing for keep their ENA files, so a partial SDL answer
    degrades to the ENA route rather than losing data.
    """
    from .ncbi import fetch_sdl, objects_to_files

    objects = fetch_sdl(session, [r.accession for r in runs])
    sample_by_run = {r.accession: r.sample for r in runs if r.sample}
    replacements = objects_to_files(objects, flavour=flavour, sample_by_run=sample_by_run)

    by_run: dict[str, list] = {}
    for f in replacements:
        owner = next(
            (o.accession for o in objects if f.url in o.links),
            None,
        )
        if owner:
            by_run.setdefault(owner, []).append(f)

    for run in runs:
        if by_run.get(run.accession):
            run.files = by_run[run.accession]
