"""Resolve raw files through NCBI's Storage Data Locator (SDL).

SDL is the service the SRA Toolkit itself queries to find data, and it is a
plain HTTP API -- so the locations it returns can be fetched with aria2
without installing the toolkit at all.

For each run it lists two kinds of object:

* ``srapub`` -- the normalised ``.sra`` archive, which ``fasterq-dump``
  converts to FASTQ.
* ``remote`` -- the file the submitter originally uploaded, often a BAM.

The ``.sra`` is frequently much smaller than the original. For SRR16036138
it is 48 GB against 119 GB for the submitted BAM, so when the goal is FASTQ
the ``.sra`` route moves roughly 2.5x less data.

SDL publishes an md5 for every object, so transfers are verifiable.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Any

import requests

from ..errors import ResolveError
from ..models import RemoteFile

SDL_ENDPOINT = "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"

#: SDL object prefixes.
SRA_OBJECT = "srapub"
ORIGINAL_OBJECT = "remote"

#: What the caller wants out of SDL.
FLAVOURS = ("sra", "original")


@dataclass
class SdlObject:
    """One object SDL knows about for a run."""

    accession: str
    name: str
    size: int | None
    md5: str | None
    object_type: str
    #: ``"sra"`` or ``"original"``.
    flavour: str
    links: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return dataclasses.asdict(self)


def _flavour(object_id: str, file_type: str) -> str:
    """Classify an SDL file record as the .sra archive or the original upload."""
    if object_id.startswith(SRA_OBJECT) or file_type.lower() == "sra":
        return "sra"
    return "original"


def parse_sdl(payload: dict[str, Any]) -> list[SdlObject]:
    """Parse an SDL ``retrieve`` response into objects.

    Runs SDL reports as failed are skipped rather than raising, so one
    restricted run does not abort a whole project.
    """
    objects: list[SdlObject] = []
    for bundle in payload.get("result") or []:
        if bundle.get("status") != 200:
            continue
        accession = bundle.get("bundle") or ""
        for record in bundle.get("files") or []:
            object_id = str(record.get("object") or "")
            links = [
                loc["link"]
                for loc in (record.get("locations") or [])
                if loc.get("link")
            ]
            if not links:
                continue
            objects.append(
                SdlObject(
                    accession=record.get("accession") or accession,
                    name=record.get("name") or links[0].rsplit("/", 1)[-1],
                    size=record.get("size"),
                    md5=record.get("md5"),
                    object_type=object_id,
                    flavour=_flavour(object_id, str(record.get("type") or "")),
                    links=links,
                )
            )
    return objects


def fetch_sdl(
    session: requests.Session, run_accessions: list[str]
) -> list[SdlObject]:
    """Ask SDL where the data for *run_accessions* lives.

    Raises:
        ResolveError: if SDL cannot be reached or returns an unusable body.
    """
    if not run_accessions:
        return []
    try:
        response = session.post(
            SDL_ENDPOINT, data={"acc": ",".join(run_accessions)}, timeout=120
        )
        response.raise_for_status()
        payload = response.json()
    except (requests.RequestException, ValueError) as exc:
        raise ResolveError(f"NCBI SDL lookup failed: {exc}") from exc
    return parse_sdl(payload)


def objects_to_files(
    objects: list[SdlObject],
    flavour: str = "sra",
    sample_by_run: dict[str, str] | None = None,
) -> list[RemoteFile]:
    """Convert SDL objects into downloadable file records.

    Args:
        flavour: ``"sra"`` for the normalised archive, ``"original"`` for the
            submitter's own upload.
        sample_by_run: optional run accession to GEO sample mapping, so files
            land under the right sample directory.
    """
    if flavour not in FLAVOURS:
        raise ValueError(f"Unknown flavour {flavour!r}. Use one of: {', '.join(FLAVOURS)}.")

    sample_by_run = sample_by_run or {}
    files: list[RemoteFile] = []
    for obj in objects:
        if obj.flavour != flavour:
            continue
        # SDL may list several mirrors; the first is the preferred location.
        url = obj.links[0]
        name = obj.name
        if flavour == "sra" and not name.lower().endswith(".sra"):
            # srapub objects are named for the run; give the file an extension
            # so downstream tools recognise it.
            name = f"{name}.sra"
        files.append(
            RemoteFile(
                url=url,
                filename=name,
                size=obj.size,
                kind="raw",
                accession=sample_by_run.get(obj.accession, obj.accession),
                md5=obj.md5,
                description=f"NCBI SDL {obj.object_type}",
            )
        )
    return files
