"""Decide where each file lands under the destination directory.

The flat layout writes every file into one directory, matching the manifest
exactly. The by-sample layout groups files under their GEO sample accession
and, for recognised 10x Cell Ranger triplets, restores the canonical
``barcodes.tsv.gz`` / ``features.tsv.gz`` / ``matrix.mtx.gz`` names -- which
is what ``scanpy.read_10x_mtx()`` and Seurat's ``Read10X()`` expect.
"""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Sequence

from .models import RemoteFile

LAYOUTS = ("flat", "by-sample")

#: A 10x triplet member, e.g. GSM5008737_RNA_3P-barcodes.tsv.gz. The GEO
#: accession and the submitter's modality label both precede the role.
TENX_RE = re.compile(
    r"^GSM\d+[_-].*?(?P<role>barcodes|features|genes|matrix)"
    r"\.(?P<ext>tsv|mtx)\.gz$",
    re.IGNORECASE,
)

#: Accessions that identify an individual sample rather than a series.
SAMPLE_RE = re.compile(r"^GSM\d+$", re.IGNORECASE)


def canonical_tenx_name(filename: str) -> str | None:
    """Return the Cell Ranger name for a 10x triplet member, else None."""
    match = TENX_RE.match(filename)
    if not match:
        return None
    return f"{match.group('role').lower()}.{match.group('ext').lower()}.gz"


def apply_layout(
    files: Sequence[RemoteFile], layout: str = "by-sample"
) -> list[RemoteFile]:
    """Set ``relpath`` on each file according to *layout*.

    Args:
        files: the resolved files, which are mutated in place and returned.
        layout: ``"by-sample"`` (the default, since per-sample directories
            are what single-cell loaders expect) or ``"flat"``.

    For ``by-sample``, files belonging to a GSM accession go into a directory
    named for that sample; everything else stays at the top level. Canonical
    10x renaming is applied only when it is unambiguous: if two files in one
    sample directory would end up with the same name -- which happens when a
    GSM holds several modalities -- the original filenames are kept for that
    sample instead, because silently overwriting data is worse than an
    inconvenient name.
    """
    if layout not in LAYOUTS:
        raise ValueError(f"Unknown layout {layout!r}. Use one of: {', '.join(LAYOUTS)}.")

    if layout == "flat":
        for f in files:
            f.relpath = f.filename
        return list(files)

    grouped: dict[str | None, list[RemoteFile]] = defaultdict(list)
    for f in files:
        accession = f.accession or ""
        grouped[accession.upper() if SAMPLE_RE.match(accession) else None].append(f)

    for sample, members in grouped.items():
        if sample is None:
            for f in members:
                f.relpath = f.filename
            continue

        proposed = {f.filename: canonical_tenx_name(f.filename) for f in members}
        renamed = [n for n in proposed.values() if n]
        collides = len(renamed) != len(set(renamed))

        for f in members:
            name = proposed[f.filename]
            f.relpath = f"{sample}/{name if name and not collides else f.filename}"

    return list(files)
