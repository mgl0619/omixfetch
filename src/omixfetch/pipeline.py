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
3. ``original``   -- the submitter's own file (often BAM) from NCBI S3, then
                     10x ``bamtofastq`` if installed.

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
    #: convert_bam(bam_path, outdir) -> list[Path]
    convert_bam: Callable[[Path, Path], list[Path]] | None = None
    #: check_bam(url) -> None, raising if bamtofastq could not convert it.
    #: Costs a few MB instead of the whole file.
    check_bam: Callable[[str], None] | None = None
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

    sample_dir = out_dir / "fastq" / (plan.sample or plan.run)
    work_dir = out_dir / "work" / plan.run
    routes = plan.available(order)
    reason = plan.incomplete_sra()
    if reason and "sra" in routes and not allow_incomplete:
        routes = [r for r in routes if r != "sra"]
        tools.log(f"  {plan.run}: skipping sra route: {reason}")
    if not routes:
        return store.update(
            plan.run, sample=plan.sample, status=SKIPPED,
            error="no public source offers this run (controlled access?)",
        )

    attempts = list(state.attempts)
    for route in routes:
        attempts.append(route)
        store.update(plan.run, sample=plan.sample, route=route,
                     expected_spots=plan.expected_spots, attempts=attempts,
                     status=PENDING, error=None)
        try:
            if check_space:
                check_disk(out_dir, disk_needed(plan, route))
            fastq, spots = _run_route(plan, route, work_dir, sample_dir, tools)
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
            spots_written=spots, error=None,
        )

    return store.get(plan.run)


def _run_route(
    plan: RunPlan, route: str, work_dir: Path, sample_dir: Path, tools: Tools
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
    if route == "original" and tools.check_bam is not None:
        bam_file = next((f for f in files if ".bam" in f.filename
                         and not f.filename.endswith(".bai")), None)
        if bam_file is not None:
            tools.check_bam(bam_file.url)

    work_dir.mkdir(parents=True, exist_ok=True)
    for f in files:
        f.relpath = f.filename
    tools.download(files, work_dir)

    if route == "sra":
        sra = next((work_dir / f.filename for f in files
                    if f.filename.endswith(".sra")), work_dir / files[0].filename)
        result = tools.convert_sra(sra, work_dir)
        produced = list(getattr(result, "files", result))
        spots = getattr(result, "spots_read", None)
    else:  # original
        bam = next((work_dir / f.filename for f in files
                    if ".bam" in f.filename and not f.filename.endswith(".bai")),
                   work_dir / files[0].filename)
        produced = tools.convert_bam(bam, work_dir)
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
        plans.append(
            RunPlan(
                run=run.accession,
                sample=run.sample,
                expected_spots=expected,
                routes={k: v for k, v in routes.items() if v},
                mates_missing=bool(ri and ri.mates_missing),
                tenx=any(o.file_type.lower() == "tenx" for o in mine),
            )
        )
    return dataset, plans
