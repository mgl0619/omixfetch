"""File selection filters applied between resolution and download."""

from __future__ import annotations

import fnmatch
from collections.abc import Sequence

from .models import RemoteFile


def apply_filters(
    files: Sequence[RemoteFile],
    *,
    include: Sequence[str] = (),
    exclude: Sequence[str] = (),
    kinds: Sequence[str] = (),
    max_size: int | None = None,
) -> list[RemoteFile]:
    """Select files by glob pattern, kind and size.

    Args:
        include: keep only files whose name matches at least one glob.
        exclude: drop files whose name matches any glob. Applied after
            ``include``, so exclusion wins.
        kinds: keep only these kinds (``supplementary``, ``metadata``, ...).
        max_size: drop files larger than this many bytes. Files with an
            unknown size are kept, since a missing size is not evidence of a
            large file.

    Returns:
        The surviving files, in their original order.
    """
    selected = list(files)

    if include:
        selected = [
            f for f in selected if any(fnmatch.fnmatch(f.filename, p) for p in include)
        ]
    if exclude:
        selected = [
            f
            for f in selected
            if not any(fnmatch.fnmatch(f.filename, p) for p in exclude)
        ]
    if kinds:
        wanted = {k.lower() for k in kinds}
        selected = [f for f in selected if f.kind.lower() in wanted]
    if max_size is not None:
        selected = [f for f in selected if f.size is None or f.size <= max_size]

    return selected


def parse_size(text: str) -> int:
    """Parse a human-readable size such as ``"500MB"`` into bytes."""
    s = text.strip().upper().replace(" ", "")
    units = {"B": 1, "K": 1 << 10, "KB": 1 << 10, "M": 1 << 20, "MB": 1 << 20,
             "G": 1 << 30, "GB": 1 << 30, "T": 1 << 40, "TB": 1 << 40}
    for suffix in sorted(units, key=len, reverse=True):
        if s.endswith(suffix):
            number = s[: -len(suffix)] or "0"
            return int(float(number) * units[suffix])
    return int(float(s))
