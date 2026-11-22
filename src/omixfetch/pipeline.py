"""Autonomous GEO -> FASTQ pipeline.

Takes a GEO series and leaves verified, gzipped FASTQ on disk, grouped by
sample. For each sequencing run it picks a route, downloads, converts if
needed, checks the result against archive metadata, and records progress so
an interrupted run resumes where it stopped.

Route preference, per run, falling through on failure:

1. ``ena-fastq``  -- ENA already serves FASTQ: download, md5-verified, done.
2. ``sra``        -- NCBI's normalised .sra from S3, then fasterq-dump.
                     Smallest transfer for aligned submissions, but needs
                     NCBI's reference service at conversion time.
3. ``original``   -- the submitter's own file (often BAM) from NCBI's AWS
                     buckets, converted back to FASTQ. The BAM header decides
                     how: ``bamtofastq`` for 10x BAMs, ``samtools fastq``
                     otherwise -- and whether the result would be complete.
                     A BAM converts back fully only if it kept the unmapped
                     reads and, for single-cell data, the barcode and UMI.

Each step is skipped when its output already exists, so rerunning the same
command is always safe.
"""

from __future__ import annotations

import json
import shutil
import threading
import time
from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from pathlib import Path

from .errors import BackendError, OmixfetchError
from .models import RemoteFile

#: Route names, in default preference order.
ROUTES = ("ena-fastq", "sra", "original")

#: Statuses a run moves through.
PENDING, DOWNLOADED, DONE, FAILED, SKIPPED = (
    "pending", "downloaded", "done", "failed", "skipped",
)

#: Free disk needed per byte transferred, by route. fasterq-dump writes
#: uncompressed FASTQ plus scratch before gzip, several times the .sra size.
#: These are deliberately generous estimates, not measurements.
DISK_FACTOR = {"ena-fastq": 1.1, "sra": 8.0, "original": 3.0}


@dataclass
class RunPlan:
    """Everything needed to fetch one run, by every available route."""

    run: str
    sample: str | None
    expected_spots: int | None
    #: Files per route; a route is only offered when it has files.
    routes: dict[str, list[RemoteFile]] = field(default_factory=dict)
    #: The library is paired but its .sra has no mates: a read was lost when
    #: NCBI loaded the submission, so converting the .sra yields half a pair.
    mates_missing: bool = False
    #: NCBI marks the original submission as a Cell Ranger (10x) BAM.
    tenx: bool = False
    #: What the submitter uploaded: "fastq", "bam", "cram" or "unknown".
    submitted: str = "unknown"
    #: BAM header assessment, for BAM/CRAM submissions (see probe.assess_bam).
    bam_verdict: object | None = None
    #: Single-cell data, judged from the BAM header or ENA library source.
    single_cell: bool = False
    #: Outcome of the decision tree; filled in before processing.
    decision: object | None = None

    def available(self, order: tuple[str, ...] = ROUTES) -> list[str]:
        return [r for r in order if self.routes.get(r)]

    def incomplete_sra(self) -> str | None:
        """Why the .sra route would produce incomplete reads, if it would.

        For a 10x BAM submission the .sra keeps only the cDNA read: cell
        barcodes and UMIs lived in BAM tags and are not loaded as sequence.
        GSE184708 shows it plainly -- declared PAIRED, 98 bases per spot,
        spots_with_mates 0.
        """
        if self.tenx:
            return ".sra of a 10x BAM lacks the barcode/UMI read (R1)"
        if self.mates_missing:
            return ".sra is declared paired-end but has no mates"
        return None


@dataclass
class RunState:
    """Persisted progress for one run."""

    run: str
    sample: str | None = None
    status: str = PENDING
    route: str | None = None
    fastq: list[str] = field(default_factory=list)
    expected_spots: int | None = None
    spots_written: int | None = None
    attempts: list[str] = field(default_factory=list)
    error: str | None = None
    updated: float = 0.0
    #: Decision-tree branch taken, and routes excluded with reasons.
    branch: str | None = None
    excluded: dict = field(default_factory=dict)
    submitted: str | None = None
    #: Wall-clock timing and bytes moved for the successful route.
    started: float | None = None
    finished: float | None = None
    bytes: int | None = None


class StateStore:
    """Thread-safe JSON record of per-run progress.

    Written atomically after every change, so a crash never leaves a
    half-written state file and a rerun sees exactly what finished.
    """

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        self._lock = threading.Lock()
        self.runs: dict[str, RunState] = {}
        if self.path.exists():
            raw = json.loads(self.path.read_text())
            self.runs = {k: RunState(**v) for k, v in raw.get("runs", {}).items()}

    def get(self, run: str) -> RunState:
        with self._lock:
            return self.runs.setdefault(run, RunState(run=run))

    def update(self, run: str, **changes) -> RunState:
        with self._lock:
            state = self.runs.setdefault(run, RunState(run=run))
            for key, value in changes.items():
                setattr(state, key, value)
            state.updated = time.time()
            self._save()
            return state

    def _save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(
            json.dumps({"runs": {k: asdict(v) for k, v in self.runs.items()}}, indent=2)
        )
        tmp.replace(self.path)


def is_complete(state: RunState) -> bool:
    """A run counts as done only if its recorded FASTQ files still exist."""
    return (
        state.status == DONE
        and bool(state.fastq)
        and all(Path(p).exists() for p in state.fastq)
    )


def disk_needed(plan: RunPlan, route: str) -> int:
    """Bytes of free space a route is expected to need for this run."""
    size = sum(f.size or 0 for f in plan.routes.get(route, []))
    return int(size * DISK_FACTOR.get(route, 2.0))


def check_disk(path: Path, needed: int) -> None:
    """Raise if *path* has less free space than *needed*."""
    path.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(path).free
    if needed and free < needed:
        raise BackendError(
            f"not enough disk: about {needed / 1024**3:.0f} GiB needed, "
            f"{free / 1024**3:.0f} GiB free at {path}"
        )


def verify_spots(expected: int | None, observed: int | None) -> str | None:
    """Return a problem description if spot counts disagree, else None.

    Unknown counts on either side are not treated as failures: a missing
    number is not evidence of a truncated file.
    """
    # Zero means the archive never counted the reads (ENA reports 0 for
    # submissions it did not process), not that the run is empty.
    if not expected or observed is None:
        return None
    if expected != observed:
        return f"spot count mismatch: expected {expected:,}, got {observed:,}"
    return None


@dataclass
class Tools:
    """External actions, injected so the pipeline is testable offline."""

    #: download(files, dest_dir) -> None
    download: Callable[[list[RemoteFile], Path], None]
    #: convert_sra(sra_path, outdir) -> ConversionResult-like
    convert_sra: Callable[[Path, Path], object] | None = None
    #: convert_bam(bam_path, outdir, converter) -> list[Path]
    convert_bam: Callable[..., list[Path]] | None = None
    #: check_bam(url) -> BamAssessment, read from the header with a range
    #: request: a few MB instead of the whole file.
    check_bam: Callable[[str], object] | None = None
    #: Converter executables actually installed, e.g. {"bamtofastq", "samtools"}.
    converters: frozenset = frozenset({"bamtofastq"})
    log: Callable[[str], None] = print


def process_run(
    plan: RunPlan,
    out_dir: Path,
    store: StateStore,
    tools: Tools,
    *,
    order: tuple[str, ...] = ROUTES,
    keep_intermediate: bool = False,
    check_space: bool = True,
    allow_incomplete: bool = False,
) -> RunState:
    """Take one run to verified FASTQ, trying routes in *order*.

    The .sra route is dropped for runs where it would yield incomplete reads
    (see :meth:`RunPlan.incomplete_sra`) unless *allow_incomplete* is set.
    """
    state = store.get(plan.run)
    if is_complete(state):
        tools.log(f"  {plan.run}: already done, skipping")
        return state

    from .decide import decide

    sample_dir = out_dir / "fastq" / (plan.sample or plan.run)
    work_dir = out_dir / "work" / plan.run
    decision = plan.decision or decide(
        plan, allow_incomplete=allow_incomplete,
        order_override=order if tuple(order) != ROUTES else None,
    )
    routes = list(decision.order)
    store.update(plan.run, branch=decision.branch, excluded=dict(decision.excluded),
                 submitted=plan.submitted)
    for route, why in decision.excluded.items():
        if plan.routes.get(route):
            tools.log(f"  {plan.run}: not using {route}: {why}")
    if not routes:
        if not plan.available():
            return store.update(
                plan.run, sample=plan.sample, status=SKIPPED,
                error="no public source offers this run (controlled access?)",
            )
        return store.update(
            plan.run, sample=plan.sample, status=FAILED,
            error=f"no faithful route: {decision.branch} (--allow-incomplete to accept loss)",
        )

    attempts = list(state.attempts)
    for route in routes:
        attempts.append(route)
        store.update(plan.run, sample=plan.sample, route=route,
                     expected_spots=plan.expected_spots, attempts=attempts,
                     status=PENDING, error=None)
        started = time.time()
        try:
            if check_space:
                check_disk(out_dir, disk_needed(plan, route))
            fastq, spots = _run_route(
                plan, route, work_dir, sample_dir, tools, allow_incomplete
            )
            problem = verify_spots(plan.expected_spots, spots)
            if problem:
                raise BackendError(problem)
        except (OmixfetchError, OSError) as exc:
            tools.log(f"  {plan.run}: route {route} failed: {str(exc)[:200]}")
            store.update(plan.run, status=FAILED, error=f"{route}: {exc}")
            continue

        if not keep_intermediate:
            shutil.rmtree(work_dir, ignore_errors=True)
        tools.log(f"  {plan.run}: done via {route} ({len(fastq)} file(s))")
        return store.update(
            plan.run, status=DONE, fastq=[str(p) for p in fastq],
            spots_written=spots, error=None, started=started, finished=time.time(),
            bytes=sum(f.size or 0 for f in plan.routes[route]),
        )

    return store.get(plan.run)


def _run_route(
    plan: RunPlan,
    route: str,
    work_dir: Path,
    sample_dir: Path,
    tools: Tools,
    allow_incomplete: bool = False,
) -> tuple[list[Path], int | None]:
    """Execute one route; return final FASTQ paths and the spot count seen."""
    files = plan.routes[route]
    sample_dir.mkdir(parents=True, exist_ok=True)

    if route == "ena-fastq":
        # FASTQ arrives finished and md5-checked by the backend.
        for f in files:
            f.relpath = f.filename
        tools.download(files, sample_dir)
        return [sample_dir / f.filename for f in files], None

    # Refuse a route whose converter is missing before downloading anything:
    # discovering it afterwards would waste a transfer of hundreds of GB.
    if route == "sra" and tools.convert_sra is None:
        raise BackendError("no fasterq-dump available for the sra route")
    if route == "original" and tools.convert_bam is None:
        raise BackendError(
            "no bamtofastq available to convert the original submission; "
            "run `omixfetch setup` to install it"
        )
    converter = "bamtofastq"
    if route == "original" and plan.bam_verdict is not None:
        verdict = plan.bam_verdict
        if verdict.converter is None:
            raise BackendError(f"cannot convert back to FASTQ: {verdict.reason}")
        if not verdict.acceptable and not allow_incomplete:
            raise BackendError(f"conversion would be {verdict.complete}: {verdict.reason}")
        converter = verdict.converter
        if converter not in tools.converters:
            raise BackendError(f"this BAM needs {converter}, which is not installed")
    elif route == "original" and tools.check_bam is not None:
        bam_file = next((f for f in files if ".bam" in f.filename
                         and not f.filename.endswith(".bai")), None)
        if bam_file is not None:
            verdict = tools.check_bam(bam_file.url)
            if verdict is not None:
                if verdict.converter is None:
                    raise BackendError(f"cannot convert back to FASTQ: {verdict.reason}")
                if not verdict.acceptable and not allow_incomplete:
                    raise BackendError(
                        f"conversion would be {verdict.complete}: {verdict.reason} "
                        f"(--allow-incomplete to accept)"
                    )
                converter = verdict.converter
        if converter not in tools.converters:
            raise BackendError(
                f"this BAM needs {converter}, which is not installed"
                + (" (conda install -c bioconda samtools)" if converter == "samtools" else "")
            )

    work_dir.mkdir(parents=True, exist_ok=True)
    for f in files:
        f.relpath = f.filename
    tools.download(files, work_dir)

    if route == "sra":
        sra = next((work_dir / f.filename for f in files
                    if f.filename.endswith(".sra")), work_dir / files[0].filename)
        result = (tools.convert_sra(sra, work_dir, single_cell=True)
                  if plan.single_cell else tools.convert_sra(sra, work_dir))
        produced = list(getattr(result, "files", result))
        spots = getattr(result, "spots_read", None)
    else:  # original
        bam = next((work_dir / f.filename for f in files
                    if ".bam" in f.filename and not f.filename.endswith(".bai")),
                   work_dir / files[0].filename)
        produced = tools.convert_bam(bam, work_dir, converter)
        spots = None

    final = []
    for path in produced:
        target = sample_dir / Path(path).name
        shutil.move(str(path), target)
        final.append(target)
    if not final:
        raise BackendError(f"{route} produced no FASTQ files")
    return final, spots


def write_manifest(store: StateStore, path: Path) -> Path:
    """Write a TSV of every run's outcome and files."""
    rows = ["\t".join(["sample", "run", "status", "route", "expected_spots",
                        "spots_written", "fastq", "error"])]
    for state in sorted(store.runs.values(), key=lambda s: (s.sample or "", s.run)):
        rows.append("\t".join([
            state.sample or "", state.run, state.status, state.route or "",
            "" if state.expected_spots is None else str(state.expected_spots),
            "" if state.spots_written is None else str(state.spots_written),
            ";".join(state.fastq), (state.error or "").replace("\t", " ")[:300],
        ]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(rows) + "\n")
    return path


#: Suffixes that mark a file as finished FASTQ rather than a submitted BAM.
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

#: SDL accepts batches; keep requests modest for large series.
SDL_BATCH = 50


def plan_runs(accession: str, session) -> tuple[object, list[RunPlan]]:
    """Resolve a GEO series into one :class:`RunPlan` per sequencing run.

    ENA supplies the run list, sample mapping and expected spot counts; NCBI
    SDL supplies the .sra and original-submission locations. A failed SDL
    lookup leaves only the ENA route rather than aborting.
    """
    from .errors import ResolveError
    from .resolve import resolve
    from .sources.ncbi import fetch_sdl, objects_to_files
    from .sources.runinfo import fetch_runinfo

    dataset = resolve(
        accession,
        session=session,
        resolver="native",
        include_samples=True,
        include_reads=True,
        fetch_sizes=False,
    )
    runs = list(dataset.runs)

    objects = []
    accessions = [r.accession for r in runs]
    for i in range(0, len(accessions), SDL_BATCH):
        try:
            objects += fetch_sdl(session, accessions[i : i + SDL_BATCH])
        except ResolveError:
            continue

    info = fetch_runinfo(session, accessions)

    plans = []
    for run in runs:
        mine = [o for o in objects if o.accession == run.accession]
        ri = info.get(run.accession)
        # NCBI's count first: ENA reports 0 for submissions it never processed.
        expected = (ri.spots if ri and ri.spots else None) or (run.read_count or None)
        routes = {
            "ena-fastq": [
                f for f in run.files if f.filename.lower().endswith(FASTQ_SUFFIXES)
            ],
            "sra": objects_to_files(mine, flavour="sra"),
            "original": objects_to_files(mine, flavour="original"),
        }
        from .decide import submission_format

        originals = [o.name for o in mine if o.flavour == "original"]
        source = (run.library_source or "").upper()
        plans.append(
            RunPlan(
                run=run.accession,
                sample=run.sample,
                expected_spots=expected,
                routes={k: v for k, v in routes.items() if v},
                mates_missing=bool(ri and ri.mates_missing),
                tenx=any(o.file_type.lower() == "tenx" for o in mine),
                submitted=submission_format(originals, run.raw.get("submitted_format")),
                single_cell="SINGLE CELL" in source,
            )
        )

    probe_bams(plans, session)
    return dataset, plans


def probe_bams(plans: list[RunPlan], session, workers: int = 8) -> None:
    """Read the header of every BAM submission (a few MB each), in parallel.

    The header is the only reliable record of whether the BAM converts back
    to complete FASTQ, and of whether the data is single-cell.
    """
    from concurrent.futures import ThreadPoolExecutor

    from .probe import IMPOSSIBLE, BamAssessment, assess_bam, bam_header

    def one(plan: RunPlan) -> None:
        bam = next((f for f in plan.routes.get("original", [])
                    if ".bam" in f.filename and not f.filename.endswith(".bai")), None)
        if bam is None:
            return
        try:
            plan.bam_verdict = assess_bam(bam_header(session, bam.url))
        except Exception as exc:  # unreadable header: recorded, decided later
            plan.bam_verdict = BamAssessment(None, "unknown", None,
                                             f"BAM header could not be read: {exc}")
            return
        if plan.bam_verdict.converter == "bamtofastq" or plan.bam_verdict.complete == IMPOSSIBLE:
            plan.single_cell = True

    targets = [p for p in plans if p.submitted in ("bam", "cram")]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(one, targets))


def measure_rates(plans: list[RunPlan], session, sample_bytes: int = 16 * 1024**2) -> dict:
    """Measure throughput once per series from each source the plan uses.

    Returns MiB/s keyed ``"ena"`` and ``"ncbi"``; a source that fails to
    respond is simply absent, and the decision falls back to a default.
    """
    from .benchmark import benchmark_source

    picks = {}
    for plan in plans:
        for route, key in (("ena-fastq", "ena"), ("sra", "ncbi"), ("original", "ncbi")):
            if key not in picks and plan.routes.get(route):
                picks[key] = max(plan.routes[route], key=lambda f: f.size or 0)
    rates = {}
    for key, f in picks.items():
        result = benchmark_source(key, f, session, sample_bytes)
        if result.mib_per_s:
            rates[key] = result.mib_per_s
    return rates


def apply_decisions(
    plans: list[RunPlan],
    rates: dict | None = None,
    allow_incomplete: bool = False,
    order_override: tuple[str, ...] | None = None,
) -> None:
    """Run the decision tree for every plan, storing the result on it."""
    from .decide import decide

    for plan in plans:
        plan.decision = decide(plan, rates, allow_incomplete, order_override)
