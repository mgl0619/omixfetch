"""Choose how to turn each sequencing run into FASTQ.

The decision is driven by what the submitter actually uploaded, because that
fixes what can be recovered:

    run
    |-- no public files ........................ skip (controlled access)
    |-- submitted as BAM/CRAM  (header read with a 4 MB range request)
    |   |-- 10x_bam_to_fastq lines ............. AWS BAM -> bamtofastq
    |   |                                         .sra and ENA FASTQ excluded: derived
    |   |                                         from the BAM, they lost the barcode read
    |   |-- single-cell, no 10x lines .......... no faithful route
    |   |-- aligner kept unmapped reads ........ AWS BAM -> samtools fastq
    |   `-- dropped unmapped / unknown ......... no faithful route
    `-- submitted as FASTQ
        `-- faithful: ENA FASTQ, AWS .sra -> fasterq-dump
            ranked by estimated time: size / measured throughput,
            plus conversion for the .sra

"No faithful route" fails the run with the reason, unless the caller accepts
incomplete output, in which case the lossy routes are tried.

Signals, in order of trust (observed on GSE184708 and GSE220167):

* NCBI SDL's original-file record: a ``.bam``/``.cram`` means a BAM/CRAM
  submission; ``TenX`` marks a Cell Ranger BAM.
* ENA's ``submitted_format``.
* The BAM header itself -- the only reliable single-cell signal. NCBI's
  LibrarySource said TRANSCRIPTOMIC for GSE184708, which is 10x data.
"""

from __future__ import annotations

from dataclasses import dataclass, field

FASTQ, BAM, CRAM, UNKNOWN_FORMAT = "fastq", "bam", "cram", "unknown"

#: Relative cost of converting one byte of .sra versus just downloading it,
#: used only to rank otherwise-faithful routes. A deliberate rough estimate:
#: fasterq-dump plus gzip roughly doubles the wall time of a fast transfer.
SRA_CONVERSION_PENALTY = 1.0

#: Throughput assumed for a source that was not measured, in MiB/s.
DEFAULT_RATE = 10.0


@dataclass
class Decision:
    """The outcome of the tree for one run."""

    #: Routes to try, best first.
    order: list[str]
    #: The branch taken, as human-readable text for the dry run and report.
    branch: str
    #: Routes deliberately not used, with the reason.
    excluded: dict[str, str] = field(default_factory=dict)
    #: Estimated hours per candidate route, when ranking was possible.
    estimates: dict[str, float] = field(default_factory=dict)
    #: True if every route left would lose data.
    lossy: bool = False


def submission_format(original_names: list[str], ena_submitted: str | None) -> str:
    """Infer what the submitter uploaded.

    NCBI lists an original file only for non-FASTQ submissions (BAM, CRAM,
    archives); FASTQ submissions are loaded directly into .sra.
    """
    names = [n.lower() for n in original_names]
    if any(n.endswith(".cram") or ".cram." in n for n in names):
        return CRAM
    if any(".bam" in n and not n.endswith(".bai") for n in names):
        return BAM
    fmt = (ena_submitted or "").upper()
    if "CRAM" in fmt:
        return CRAM
    if "BAM" in fmt:
        return BAM
    if "FASTQ" in fmt or not names:
        return FASTQ
    return UNKNOWN_FORMAT


def _hours(size: int, rate_mib: float) -> float:
    return size / (max(rate_mib, 0.01) * 1024**2) / 3600


def decide(
    plan,
    rates: dict[str, float] | None = None,
    allow_incomplete: bool = False,
    order_override: tuple[str, ...] | None = None,
) -> Decision:
    """Walk the tree for one :class:`~omixfetch.pipeline.RunPlan`.

    Args:
        rates: measured MiB/s per source, keys ``"ena"`` and ``"ncbi"``.
        allow_incomplete: let lossy routes run after faithful ones fail.
        order_override: an explicit route order from the user. Exclusions
            for data loss still apply unless *allow_incomplete* is set.
    """
    rates = rates or {}
    available = [r for r in ("ena-fastq", "sra", "original") if plan.routes.get(r)]
    if not available:
        return Decision([], "no public files (controlled access?)")

    excluded: dict[str, str] = {}
    verdict = plan.bam_verdict
    fmt = plan.submitted

    if fmt in (BAM, CRAM):
        derived = "derived from the submitted BAM"
        if verdict is not None and getattr(verdict, "converter", None) == "bamtofastq":
            branch = "BAM submission, 10x: AWS BAM -> bamtofastq"
            for r in ("sra", "ena-fastq"):
                excluded[r] = f"{derived}; lost the barcode/UMI read (R1)"
            faithful = ["original"]
        elif verdict is not None and getattr(verdict, "acceptable", False):
            branch = f"BAM submission, {verdict.aligner} kept unmapped reads: AWS BAM -> samtools"
            faithful = ["original"]
            for r in ("sra", "ena-fastq"):
                if plan.mates_missing:
                    excluded[r] = f"{derived}; paired library with no mates"
            faithful += [r for r in ("sra", "ena-fastq") if r not in excluded]
        else:
            reason = (verdict.reason if verdict is not None
                      else "BAM header could not be read")
            branch = f"BAM submission, no faithful route: {reason}"
            excluded["original"] = reason
            for r in ("sra", "ena-fastq"):
                excluded[r] = f"{derived}; carries the same loss"
            faithful = []
    else:
        branch = "FASTQ submission: fastest faithful route"
        faithful = ["ena-fastq", "sra", "original"]
        if plan.mates_missing:
            excluded["sra"] = ".sra declared paired-end but has no mates"
            excluded["ena-fastq"] = "generated from the same incomplete .sra"

    faithful = [r for r in faithful if r in available and r not in excluded]

    estimates: dict[str, float] = {}
    for route in faithful:
        size = sum(f.size or 0 for f in plan.routes[route])
        source = "ena" if route == "ena-fastq" else "ncbi"
        hours = _hours(size, rates.get(source, DEFAULT_RATE))
        needs_conversion = route == "sra" or (
            route == "original"
            and not all(f.filename.lower().endswith((".fastq.gz", ".fq.gz"))
                        for f in plan.routes[route])
        )
        if needs_conversion:
            hours *= 1 + SRA_CONVERSION_PENALTY
        estimates[route] = hours
    # Within a branch, the first faithful route is a no-conversion route only
    # when it is also fastest; ranking by estimate makes that explicit.
    if fmt not in (BAM, CRAM):
        # Stable sort: equal estimates keep the tree's preference order.
        faithful.sort(key=lambda r: estimates.get(r, float("inf")))

    if order_override:
        faithful = [r for r in order_override if r in faithful]

    lossy_tail = []
    if allow_incomplete:
        lossy_tail = [r for r in (order_override or ("original", "sra", "ena-fastq"))
                      if r in available and r in excluded and r not in faithful]

    return Decision(
        order=faithful + lossy_tail,
        branch=branch,
        excluded=excluded,
        estimates=estimates,
        lossy=not faithful and bool(lossy_tail),
    )
