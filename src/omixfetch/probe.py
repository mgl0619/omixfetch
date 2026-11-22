"""Cheap checks on remote files, made before committing to a download."""

from __future__ import annotations

import zlib

import requests

#: Bytes fetched to read a BAM header. Large reference dictionaries can push
#: the header past 1 MB; 4 MB has been ample in practice.
HEADER_BYTES = 4 * 1024 * 1024

TENX_MARKER = "10x_bam_to_fastq:"


def bam_header(session: requests.Session, url: str, nbytes: int = HEADER_BYTES) -> str:
    """Return the SAM text header of a remote BAM, fetched by range request.

    BAM is BGZF: a series of gzip members. Only the first few are inflated --
    enough to recover the header -- so this costs megabytes, not the file.

    Raises:
        ValueError: if the bytes are not a BAM.
    """
    raw = session.get(url, headers={"Range": f"bytes=0-{nbytes - 1}"}, timeout=120).content
    out, pos = b"", 0
    while pos < len(raw):
        inflater = zlib.decompressobj(31)
        try:
            out += inflater.decompress(raw[pos:])
        except zlib.error:
            break  # a truncated final member; the header is already out
        used = len(raw) - pos - len(inflater.unused_data)
        if used <= 0:
            break
        pos += used
        if len(out) >= 8 and len(out) >= 8 + int.from_bytes(out[4:8], "little"):
            break
    if out[:4] != b"BAM\x01":
        raise ValueError("not a BAM file")
    length = int.from_bytes(out[4:8], "little")
    return out[8 : 8 + length].decode("utf-8", "replace")


def is_tenx_bam(header: str) -> bool:
    """True if the header carries the @CO lines 10x bamtofastq needs."""
    return any(
        line.startswith("@CO") and TENX_MARKER in line for line in header.splitlines()
    )
