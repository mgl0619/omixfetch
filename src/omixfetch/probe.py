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


# ---------------------------------------------------------------------------
# Can this BAM be turned back into the FASTQ that was sequenced?
# ---------------------------------------------------------------------------

import re  # noqa: E402
import shlex  # noqa: E402
from dataclasses import dataclass  # noqa: E402

#: Values of BamAssessment.complete.
COMPLETE, MAPPED_ONLY, UNKNOWN, IMPOSSIBLE = "complete", "mapped-only", "unknown", "impossible"

_SINGLE_CELL_HINTS = ("cellranger", "sc_rna_counter", "starsolo", "--soloType",
                      "dropseq", "alevin", "10x", "CB:Z", "UB:Z")


@dataclass
class BamAssessment:
    """Whether and how a BAM can be converted back to FASTQ."""

    #: ``"bamtofastq"``, ``"samtools"`` or None when no faithful conversion exists.
    converter: str | None
    #: One of COMPLETE, MAPPED_ONLY, UNKNOWN, IMPOSSIBLE.
    complete: str
    #: The aligner that wrote the BAM, as recorded in @PG, if known.
    aligner: str | None
    reason: str

    @property
    def acceptable(self) -> bool:
        """Safe to convert without the user explicitly accepting a loss."""
        return self.converter is not None and self.complete == COMPLETE


def _programs(header: str) -> list[tuple[str, str]]:
    """Return ``(program name, command line)`` for each @PG line."""
    out = []
    for line in header.splitlines():
        if not line.startswith("@PG"):
            continue
        fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
        out.append(((fields.get("PN") or fields.get("ID") or "").lower(), fields.get("CL", "")))
    return out


def _star_keeps_unmapped(cl: str) -> bool:
    try:
        tokens = shlex.split(cl)
    except ValueError:
        tokens = cl.split()
    for i, token in enumerate(tokens):
        if token == "--outSAMunmapped":
            return i + 1 < len(tokens) and tokens[i + 1] == "Within"
    return False  # STAR's default is None: unmapped reads are not written


def assess_bam(header: str) -> BamAssessment:
    """Decide from a BAM's header whether it converts back to complete FASTQ.

    The rule: a BAM converts back fully only if it kept the unmapped reads
    and, for single-cell data, the barcode and UMI. 10x writes explicit
    ``@CO 10x_bam_to_fastq`` instructions for the latter; otherwise the
    aligner's recorded command line says whether unmapped reads were kept.
    """
    if is_tenx_bam(header):
        return BamAssessment("bamtofastq", COMPLETE, "cellranger",
                             "10x BAM: bamtofastq rebuilds R1/R2/I1 from its tags")

    programs = _programs(header)
    text = header.lower()
    single_cell = any(h.lower() in text for h in _SINGLE_CELL_HINTS)
    if single_cell:
        return BamAssessment(
            None, IMPOSSIBLE, programs[0][0] if programs else None,
            "single-cell BAM without 10x_bam_to_fastq header lines: cell "
            "barcodes and UMIs cannot be rebuilt as reads",
        )

    for name, cl in programs:
        if name == "star":
            if _star_keeps_unmapped(cl):
                return BamAssessment("samtools", COMPLETE, "star",
                                     "STAR with --outSAMunmapped Within keeps every read")
            return BamAssessment("samtools", MAPPED_ONLY, "star",
                                 "STAR without --outSAMunmapped Within drops unmapped reads")
        if name in ("hisat2", "bowtie2"):
            if re.search(r"(^|\s)--no-unal(\s|$)", cl):
                return BamAssessment("samtools", MAPPED_ONLY, name,
                                     f"{name} ran with --no-unal, dropping unmapped reads")
            return BamAssessment("samtools", COMPLETE, name,
                                 f"{name} keeps unmapped reads by default")
        if name in ("bwa", "minimap2"):
            if name == "minimap2" and "--sam-hit-only" in cl:
                return BamAssessment("samtools", MAPPED_ONLY, name,
                                     "minimap2 ran with --sam-hit-only")
            return BamAssessment("samtools", COMPLETE, name,
                                 f"{name} keeps unmapped reads by default")

    return BamAssessment(
        "samtools", UNKNOWN, programs[0][0] if programs else None,
        "aligner not recognised from @PG; whether unmapped reads were kept is unknown",
    )
