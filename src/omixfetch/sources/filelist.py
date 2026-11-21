"""Parse GEO's ``filelist.txt`` manifest.

Every series supplementary directory carries a ``filelist.txt`` giving the
name, size and type of each file, including the members packed inside
``<GSE>_RAW.tar``. It is worth one request: it supplies file sizes without a
HEAD per file, and it reveals that the tar is a bundle of the per-sample
supplementary files rather than additional data.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass

import requests

from ..http import get_text

#: Rows describing the archive itself use this marker in the first column.
ARCHIVE_MARKER = "Archive"


@dataclass
class FileListEntry:
    """One row of ``filelist.txt``."""

    name: str
    size: int | None
    type: str | None
    #: True for the ``_RAW.tar`` row itself rather than a member file.
    is_archive: bool = False


def parse_filelist(text: str) -> list[FileListEntry]:
    """Parse ``filelist.txt`` into entries.

    The file is tab separated with a ``#``-prefixed header
    (``#Archive/File``, ``Name``, ``Time``, ``Size``, ``Type``).
    """
    stripped = text.lstrip()
    if stripped.startswith("#"):
        stripped = stripped[1:]
    rows = list(csv.DictReader(io.StringIO(stripped), delimiter="\t"))

    entries: list[FileListEntry] = []
    for row in rows:
        name = (row.get("Name") or "").strip()
        if not name:
            continue
        raw_size = (row.get("Size") or "").strip()
        entries.append(
            FileListEntry(
                name=name,
                size=int(raw_size) if raw_size.isdigit() else None,
                type=(row.get("Type") or "").strip() or None,
                is_archive=(row.get("Archive/File") or "").strip() == ARCHIVE_MARKER,
            )
        )
    return entries


def fetch_filelist(
    session: requests.Session, suppl_base_url: str
) -> list[FileListEntry]:
    """Fetch and parse ``filelist.txt`` for a series supplementary directory.

    Returns an empty list when the series publishes no manifest, which is
    common for older or array-based submissions.
    """
    url = suppl_base_url.rstrip("/") + "/filelist.txt"
    try:
        text = get_text(session, url, expect_text=False, tries=2)
    except Exception:
        return []
    if not text.strip():
        return []
    try:
        return parse_filelist(text)
    except (csv.Error, ValueError):
        return []


def archive_members(entries: list[FileListEntry]) -> dict[str, set[str]]:
    """Return ``{archive_name: {member filenames}}``.

    ``filelist.txt`` does not nest members under their archive, so every
    non-archive row is treated as a member of the single archive the manifest
    describes. Series with no archive row yield an empty mapping.
    """
    archives = [e.name for e in entries if e.is_archive]
    if not archives:
        return {}
    members = {e.name for e in entries if not e.is_archive}
    return {archive: set(members) for archive in archives}
