"""Sequencing projects addressed directly by read-archive accession.

Many datasets never pass through GEO: a BioProject such as PRJNA672658 has
runs in SRA/ENA but no GSE. This resolver builds the same
:class:`~omixfetch.models.Dataset` the GEO resolver does, from ENA alone:

* runs and their files come from the ENA portal ``filereport`` -- the same
  call the GEO route uses, so md5s, Aspera paths and ``--source ncbi`` behave
  identically;
* study metadata comes from the portal's ``study`` search;
* samples are the BioSamples the runs belong to, with characteristics from
  the EBI BioSamples service, which mirrors NCBI's.

Runs are grouped under their BioSample accession, which is unique and stable
where submitter aliases often are not. ``samples.tsv`` maps each accession to
its alias and title.
"""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import requests

from ..errors import ResolveError
from ..http import get_json
from ..models import Dataset, DatasetMetadata, SampleMetadata
from .ena import RunMetadata, _use_ncbi_locations, fetch_runs

SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
BIOSAMPLES = "https://www.ebi.ac.uk/biosamples/samples/{}"

STUDY_FIELDS = [
    "study_accession",
    "secondary_study_accession",
    "study_title",
    "study_description",
    "center_name",
    "first_public",
    "last_updated",
]

#: BioSamples attributes that restate INSDC bookkeeping rather than biology.
_BOOKKEEPING = {
    "External Id",
    "INSDC center name",
    "INSDC first public",
    "INSDC last update",
    "INSDC secondary accession",
    "INSDC status",
    "SRA accession",
    "NCBI submission model",
    "NCBI submission package",
}


def fetch_study(session: requests.Session, accession: str) -> dict[str, str]:
    """Return ENA's study record for a BioProject or SRA study accession.

    Returns an empty dict when ENA has no study under that accession, which
    is not an error by itself: runs may still be found.
    """
    column = (
        "study_accession" if accession.upper().startswith("PRJ") else "secondary_study_accession"
    )
    records = get_json(
        session,
        SEARCH,
        params={
            "result": "study",
            "query": f'{column}="{accession}"',
            "fields": ",".join(STUDY_FIELDS),
            "format": "json",
        },
    )
    return records[0] if isinstance(records, list) and records else {}


def fetch_characteristics(session: requests.Session, biosample: str) -> dict[str, str]:
    """Return a BioSample's biological attributes, or ``{}`` if unavailable.

    Attributes with several values are joined with ``"; "``. A failed lookup
    returns an empty dict: characteristics enrich the metadata but are not
    needed to download anything.
    """
    try:
        record = get_json(
            session, BIOSAMPLES.format(biosample), headers={"Accept": "application/json"}
        )
    except ResolveError:
        return {}
    out: dict[str, str] = {}
    for key, values in (record.get("characteristics") or {}).items():
        if key in _BOOKKEEPING:
            continue
        texts = [str(v.get("text")) for v in values if v.get("text") not in (None, "")]
        if texts:
            out[key] = "; ".join(texts)
    return out


def _samples_from_runs(runs: list[RunMetadata]) -> list[SampleMetadata]:
    """One sample per BioSample, in the order runs first mention them."""
    samples: dict[str, SampleMetadata] = {}
    for run in runs:
        key = run.sample_accession or run.secondary_sample_accession
        if not key:
            continue
        sample = samples.get(key)
        if sample is None:
            raw = run.raw
            sample = SampleMetadata(
                accession=key,
                title=raw.get("sample_title") or raw.get("sample_alias"),
                organism=raw.get("scientific_name"),
                source_name=raw.get("sample_alias"),
                library_strategy=run.library_strategy,
                library_source=run.library_source,
                library_selection=run.library_selection,
                instrument_model=run.instrument_model,
                relations={"BioSample": key},
                raw={
                    k: v
                    for k, v in raw.items()
                    if k in ("sample_alias", "sample_title", "scientific_name", "tax_id")
                },
            )
            if run.secondary_sample_accession:
                sample.relations["SRA sample"] = run.secondary_sample_accession
            samples[key] = sample
        sample.relations["runs"] = ";".join(
            filter(None, [sample.relations.get("runs"), run.accession])
        )
        run.sample = key
    return list(samples.values())


def resolve_project(
    accession: str,
    session: requests.Session,
    *,
    include_samples: bool = True,
    prefer: str = "fastq",
    source: str = "ena",
    sra_flavour: str = "sra",
    sample_workers: int = 8,
    **_ignored: Any,
) -> Dataset:
    """Resolve a read-archive accession into a :class:`Dataset` of raw reads.

    Args:
        accession: a BioProject, study, experiment, run or BioSample
            accession (normalised upper case).
        include_samples: fetch per-sample characteristics from BioSamples
            (one small request per sample, run concurrently).
        prefer, source, sra_flavour: as for the GEO route; see
            :func:`omixfetch.sources.ena.attach_reads`.

    GEO-specific options (series matrix, ``_RAW.tar`` handling, HEAD sizing)
    do not apply and are accepted and ignored, so callers can pass one set of
    options to any resolver. Reads are always resolved: they are the whole
    dataset here.

    Returns a dataset with no runs when ENA serves none -- usually controlled
    access -- rather than raising, so the caller can explain why.
    """
    runs = fetch_runs(session, accession, prefer=prefer)

    project = next((r.study_accession for r in runs if r.study_accession), None)
    study_key = project or (
        accession if accession.startswith(("PRJ", "SRP", "ERP", "DRP")) else None
    )
    try:
        study = fetch_study(session, study_key) if study_key else {}
    except ResolveError:
        study = {}

    samples = _samples_from_runs(runs)
    if include_samples and samples:
        with ThreadPoolExecutor(max_workers=max(1, sample_workers)) as pool:
            found = pool.map(lambda s: fetch_characteristics(session, s.accession), samples)
            for sample, chars in zip(samples, found, strict=True):
                sample.characteristics = chars

    if source == "ncbi" and runs:
        _use_ncbi_locations(runs, session, sra_flavour)

    by_sample = {s.accession: s for s in samples}
    for run in runs:
        for f in run.files:
            f.accession = run.sample or run.accession
        if run.sample in by_sample:
            by_sample[run.sample].files.extend(run.files)

    relations: dict[str, str] = {}
    bioproject = study.get("study_accession") or project
    srp = study.get("secondary_study_accession") or next(
        (
            r.raw.get("secondary_study_accession")
            for r in runs
            if r.raw.get("secondary_study_accession")
        ),
        None,
    )
    if bioproject:
        relations["BioProject"] = bioproject
    if srp:
        relations["SRA"] = srp

    organisms = sorted({s.organism for s in samples if s.organism})
    platforms = sorted({r.instrument_model for r in runs if r.instrument_model})
    metadata = DatasetMetadata(
        accession=accession,
        archive="sra",
        title=study.get("study_title") or None,
        summary=study.get("study_description") or None,
        organisms=organisms,
        platforms=platforms,
        contributors=[study["center_name"]] if study.get("center_name") else [],
        release_date=study.get("first_public") or None,
        last_update_date=study.get("last_updated") or None,
        relations=relations,
        n_samples=len(samples),
        raw=dict(study),
    )

    dataset = Dataset(metadata=metadata, samples=samples)
    dataset.runs = runs
    dataset.files = [f for run in runs for f in run.files]
    dataset.raw_documents = {
        f"{accession}.ena_study.json": json.dumps(study, indent=2) + "\n",
        f"{accession}.ena_runs.tsv": _runs_tsv(runs),
    }
    return dataset


def _runs_tsv(runs: list[RunMetadata]) -> str:
    """ENA's filereport rows as returned, so the archive record is kept."""
    if not runs:
        return ""
    columns: list[str] = []
    for run in runs:
        for key in run.raw:
            if key not in columns:
                columns.append(key)
    lines = ["\t".join(columns)]
    lines += ["\t".join(run.raw.get(c, "") for c in columns) for run in runs]
    return "\n".join(lines) + "\n"
