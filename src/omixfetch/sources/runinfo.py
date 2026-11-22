"""Per-run read layout from NCBI's SRA RunInfo.

ENA's spot counts are zero for submissions it never processed into reads,
which is common for BAM submissions; NCBI's are populated. RunInfo also
records whether spots have mates -- the signal that a paired library lost a
read on loading, as 10x BAMs do.
"""

from __future__ import annotations

import csv
import io
import time
from dataclasses import dataclass

import requests

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH = 100


@dataclass
class RunInfo:
    run: str
    spots: int | None
    bases: int | None
    layout: str | None
    spots_with_mates: int | None
    avg_length: int | None

    @property
    def mates_missing(self) -> bool:
        """Declared paired-end, yet no spot has a mate."""
        return (self.layout or "").upper() == "PAIRED" and self.spots_with_mates == 0


def _int(value: str | None) -> int | None:
    value = (value or "").strip()
    return int(value) if value.isdigit() else None


def parse_runinfo(text: str) -> dict[str, RunInfo]:
    info = {}
    for row in csv.DictReader(io.StringIO(text)):
        run = (row.get("Run") or "").strip()
        if not run or run == "Run":
            continue
        info[run] = RunInfo(
            run=run,
            spots=_int(row.get("spots")),
            bases=_int(row.get("bases")),
            layout=(row.get("LibraryLayout") or "").strip() or None,
            spots_with_mates=_int(row.get("spots_with_mates")),
            avg_length=_int(row.get("avgLength")),
        )
    return info


def fetch_runinfo(session: requests.Session, runs: list[str]) -> dict[str, RunInfo]:
    """Fetch RunInfo for *runs*; missing or failed batches are simply absent."""
    info: dict[str, RunInfo] = {}
    for i in range(0, len(runs), BATCH):
        batch = ",".join(runs[i : i + BATCH])
        for attempt in range(4):
            try:
                response = session.get(
                    EFETCH,
                    params={"db": "sra", "id": batch, "rettype": "runinfo", "retmode": "text"},
                    timeout=120,
                )
            except requests.RequestException:
                response = None
            if response is not None and response.ok and response.text.lstrip().startswith("Run"):
                info.update(parse_runinfo(response.text))
                break
            time.sleep(2 * (attempt + 1))
    return info
