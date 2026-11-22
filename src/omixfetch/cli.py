"""Command-line interface for omixfetch."""

from __future__ import annotations

import subprocess
import time
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from . import __version__
from .backends import (  # noqa: I001
    Aria2Backend,
    AsperaBackend,
    RequestsBackend,
    S5cmdBackend,
    aria2_available,
    aria2_version,
    ascp_available,
    s3_uri,
    s5cmd_available,
    s5cmd_version,
)
from .benchmark import benchmark_source, format_table
from .convert import (
    DEFAULT_CONVERT_ENV,
    ENV_NAME_VAR,
    PATH_ENV_VAR,
    create_convert_env,
    discover,
    fasterq_dump_available,
    fasterq_dump_version,
    sra_to_fastq,
)
from .errors import BackendError, OmixfetchError
from .filters import apply_filters, parse_size
from .http import make_session
from .layout import LAYOUTS, apply_layout
from .manifest import human_bytes, total_size, write_metadata
from .models import Dataset
from .pep import write_pep
from .resolve import resolve

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Download public genomics datasets from NCBI GEO in parallel.",
)
console = Console()
error_console = Console(stderr=True)

ACCESSIONS = typer.Argument(..., help="GEO accessions (GSE, GSM or GPL).")
RESOLVER = typer.Option(
    "auto",
    "--resolver",
    help="Metadata resolver: 'auto' (geofetch for series when installed), "
    "'native' or 'geofetch'.",
)
LAYOUT = typer.Option(
    "by-sample",
    "--layout",
    help="File layout: 'by-sample' (default; a directory per GSM, with 10x "
    "triplets renamed for scanpy/Seurat) or 'flat' (one directory, names as "
    "published).",
)
VERBOSE = typer.Option(
    False,
    "--verbose",
    help="Show the resolver's own progress logging (geofetch is quiet by default).",
)
SOURCE = typer.Option(
    "ena",
    "--source",
    help="Where raw-read bytes come from: 'ena' (EBI) or 'ncbi' (NCBI SDL / "
    "S3 mirror, often faster and can serve the smaller .sra archive).",
)
SRA_FLAVOUR = typer.Option(
    "sra",
    "--sra-flavour",
    help="With --source ncbi: 'sra' for the normalised archive (smaller, "
    "convert with fasterq-dump) or 'original' for the submitter's upload.",
)
CONVERT = typer.Option(
    None,
    "--convert",
    help="After downloading, convert .sra archives to this format. Currently "
    "'fastq' (via fasterq-dump from the SRA Toolkit).",
)
KEEP_SRA = typer.Option(
    False,
    "--keep-sra",
    help="Keep the .sra archives after converting. They are removed by default.",
)
THREADS = typer.Option(
    4, "--threads", "-t", help="Threads for conversion."
)
CONVERT_ENV = typer.Option(
    None,
    "--convert-env",
    help="Conda environment holding sra-tools. Found automatically if "
    f"'omixfetch setup' built it; also read from ${ENV_NAME_VAR}.",
)
FASTERQ_DUMP_PATH = typer.Option(
    None,
    "--fasterq-dump",
    help="Path to a fasterq-dump binary (or the directory holding it). It "
    f"need not be in omixfetch's environment. Also read from ${PATH_ENV_VAR}.",
)
PREFER = typer.Option(
    "fastq",
    "--prefer",
    help="For raw reads: 'fastq' (ENA-generated) or 'submitted' (the "
    "submitter's own files, often BAM/CRAM).",
)


def _resolve_all(accessions: list[str], **kwargs) -> list[Dataset]:
    """Resolve every accession, reporting failures without aborting the rest."""
    session = make_session()
    datasets: list[Dataset] = []
    for accession in accessions:
        try:
            with console.status(f"Resolving {accession} ..."):
                datasets.append(resolve(accession, session=session, **kwargs))
        except OmixfetchError as exc:
            error_console.print(f"[red]error[/red] {accession}: {exc}")
    if not datasets:
        raise typer.Exit(code=1)
    return datasets


def _version_callback(value: bool) -> None:
    if value:
        console.print(f"omixfetch {__version__}")
        console.print(f"aria2c: {aria2_version() or 'not installed'}")
        console.print(f"s5cmd:  {s5cmd_version() or 'not installed'}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        False, "--version", callback=_version_callback, is_eager=True,
        help="Show the omixfetch and aria2c versions and exit.",
    ),
) -> None:
    """omixfetch -- parallel, resumable access to public genomics datasets."""


@app.command("info")
def info(
    accessions: list[str] = ACCESSIONS,
    resolver: str = RESOLVER,
    verbose: bool = VERBOSE,
) -> None:
    """Print study-level metadata for each accession."""
    for dataset in _resolve_all(
        accessions,
        resolver=resolver,
        verbose=verbose,
        include_samples=True,
        fetch_sizes=False,
    ):
        m = dataset.metadata
        table = Table(title=f"{m.accession}  ({m.archive})", show_header=False, box=None)
        table.add_column(style="bold cyan", no_wrap=True)
        table.add_column(overflow="fold")
        for label, value in [
            ("Title", m.title),
            ("Organisms", ", ".join(m.organisms)),
            ("Study type", ", ".join(m.study_type)),
            ("Platforms", ", ".join(m.platforms)),
            ("Samples", str(m.n_samples)),
            ("Files", str(len(dataset.files))),
            ("Runs (ENA)", str(len(dataset.runs)) if dataset.runs else ""),
            ("Resolver", m.resolved_by),
            ("Total size", human_bytes(total_size(dataset.files))),
            ("Released", m.release_date),
            ("PubMed", ", ".join(m.pubmed_ids)),
            ("Contributors", ", ".join(m.contributors[:6])),
            ("Summary", (m.summary or "")[:600]),
        ]:
            if value:
                table.add_row(label, value)
        console.print(table)
        console.print()


@app.command("list")
def list_files(
    accessions: list[str] = ACCESSIONS,
    include: list[str] = typer.Option(
        [], "--include", "-i", help="Keep files matching this glob."
    ),
    exclude: list[str] = typer.Option(
        [], "--exclude", "-e", help="Drop files matching this glob."
    ),
    kind: list[str] = typer.Option(
        [], "--kind", help="Keep only these file kinds."
    ),
    max_size: str | None = typer.Option(
        None, "--max-size", help="Drop files larger than this, e.g. 500MB."
    ),
    samples: bool = typer.Option(
        True, "--samples/--no-samples", help="Include per-sample files."
    ),
    sizes: bool = typer.Option(
        True, "--sizes/--no-sizes", help="Query file sizes (GEO needs one HEAD per file)."
    ),
    reads: bool = typer.Option(
        False, "--reads", help="Also resolve raw sequencing reads via ENA."
    ),
    prefer: str = PREFER,
    source: str = SOURCE,
    sra_flavour: str = SRA_FLAVOUR,
    include_archives: bool = typer.Option(
        False,
        "--include-archives",
        help="Keep _RAW.tar even when its member files are listed separately.",
    ),
    resolver: str = RESOLVER,
    verbose: bool = VERBOSE,
) -> None:
    """List the files an accession resolves to, without downloading."""
    limit = parse_size(max_size) if max_size else None
    for dataset in _resolve_all(
        accessions,
        resolver=resolver,
        verbose=verbose,
        include_samples=samples,
        fetch_sizes=sizes,
        include_reads=reads,
        include_archives=include_archives,
        prefer=prefer,
        source=source,
        sra_flavour=sra_flavour,
    ):
        selected = apply_filters(
            dataset.files, include=include, exclude=exclude, kinds=kind, max_size=limit
        )
        table = Table(title=f"{dataset.accession}: {len(selected)} file(s)")
        table.add_column("Accession", style="dim")
        table.add_column("File")
        table.add_column("Kind", style="magenta")
        table.add_column("Size", justify="right")
        for f in selected:
            table.add_row(f.accession or "", f.filename, f.kind, human_bytes(f.size))
        console.print(table)
        console.print(f"[bold]Total:[/bold] {human_bytes(total_size(selected))}\n")


@app.command("metadata")
def metadata_only(
    accessions: list[str] = ACCESSIONS,
    outdir: Path = typer.Option(Path("."), "--outdir", "-o", help="Destination directory."),
    samples: bool = typer.Option(
        True, "--samples/--no-samples", help="Include per-sample annotation."
    ),
    reads: bool = typer.Option(
        False, "--reads", help="Also resolve ENA run metadata."
    ),
    pep: bool = typer.Option(
        False, "--pep", help="Additionally export PEP (config YAML + sample CSV)."
    ),
    resolver: str = RESOLVER,
    verbose: bool = VERBOSE,
) -> None:
    """Fetch and write metadata only, downloading no data files."""
    for dataset in _resolve_all(
        accessions,
        resolver=resolver,
        verbose=verbose,
        include_samples=samples,
        fetch_sizes=False,
        include_reads=reads,
    ):
        directory = outdir / dataset.accession / "metadata"
        written = write_metadata(dataset, directory)
        if pep:
            written.update(write_pep(dataset, directory / "pep"))
        console.print(
            f"[green]{dataset.accession}[/green]: "
            f"{dataset.metadata.n_samples} sample(s), "
            f"{len(dataset.files)} file record(s) -> {directory}"
        )
        for name, path in written.items():
            console.print(f"    {name:<8} {path.name}")


@app.command("manifest")
def manifest(
    accessions: list[str] = ACCESSIONS,
    out: Path = typer.Option(..., "--out", "-o", help="Path to write the manifest to."),
    fmt: str = typer.Option("aria2", "--format", "-f", help="aria2, tsv or json."),
    outdir: Path = typer.Option(
        Path("."), "--outdir", help="Directory the aria2 manifest should target."
    ),
    include: list[str] = typer.Option([], "--include", "-i"),
    exclude: list[str] = typer.Option([], "--exclude", "-e"),
    max_size: str | None = typer.Option(None, "--max-size"),
    resolver: str = RESOLVER,
) -> None:
    """Write a download manifest without downloading.

    The aria2 format can be fed straight to ``aria2c -i``, so the resolution
    step stays useful even inside a pipeline that manages its own transfers.
    """
    import json as _json

    from .backends.aria2 import render_input_file
    from .manifest import write_files_table

    limit = parse_size(max_size) if max_size else None
    all_files = []
    for dataset in _resolve_all(accessions, resolver=resolver, include_samples=True):
        all_files += [
            (dataset, f)
            for f in apply_filters(
                dataset.files, include=include, exclude=exclude, max_size=limit
            )
        ]

    out.parent.mkdir(parents=True, exist_ok=True)
    if fmt == "aria2":
        chunks = [
            render_input_file([f], (outdir / d.accession).resolve())
            for d, f in all_files
        ]
        out.write_text("".join(chunks))
    elif fmt == "tsv":
        write_files_table([f for _, f in all_files], out)
    elif fmt == "json":
        out.write_text(_json.dumps([f.to_dict() for _, f in all_files], indent=2))
    else:
        error_console.print(f"[red]Unknown format {fmt!r}. Use aria2, tsv or json.[/red]")
        raise typer.Exit(code=2)

    console.print(f"Wrote {len(all_files)} entries to {out} ({fmt}).")


@app.command("get")
def get(
    accessions: list[str] = ACCESSIONS,
    outdir: Path = typer.Option(Path("."), "--outdir", "-o", help="Destination directory."),
    jobs: int = typer.Option(4, "--jobs", "-j", help="Files downloaded concurrently."),
    connections: int = typer.Option(
        8, "--connections", "-x", help="Connections per file."
    ),
    backend: str = typer.Option(
        "auto",
        "--backend",
        "-b",
        help="auto, aria2, s5cmd, aspera or python. 'auto' prefers s5cmd for "
        "NCBI S3 sources when installed, then aria2.",
    ),
    include: list[str] = typer.Option(
        [], "--include", "-i", help="Keep files matching this glob."
    ),
    exclude: list[str] = typer.Option(
        [], "--exclude", "-e", help="Drop files matching this glob."
    ),
    kind: list[str] = typer.Option(
        [], "--kind", help="Keep only these file kinds."
    ),
    max_size: str | None = typer.Option(None, "--max-size", help="Skip files larger than this."),
    samples: bool = typer.Option(
        True, "--samples/--no-samples", help="Include per-sample files."
    ),
    metadata: bool = typer.Option(
        True, "--metadata/--no-metadata", help="Write metadata alongside the data."
    ),
    matrix: bool = typer.Option(
        False, "--matrix", help="Also fetch the GEO series matrix file."
    ),
    reads: bool = typer.Option(
        False, "--reads", help="Also fetch raw sequencing reads (FASTQ) via ENA."
    ),
    reads_only: bool = typer.Option(
        False, "--reads-only", help="Fetch only raw reads, no processed files."
    ),
    prefer: str = PREFER,
    source: str = SOURCE,
    sra_flavour: str = SRA_FLAVOUR,
    convert: str | None = CONVERT,
    keep_sra: bool = KEEP_SRA,
    threads: int = THREADS,
    fasterq_dump: str | None = FASTERQ_DUMP_PATH,
    convert_env: str | None = CONVERT_ENV,
    include_archives: bool = typer.Option(
        False,
        "--include-archives",
        help="Keep _RAW.tar even when its member files are listed separately.",
    ),
    pep: bool = typer.Option(
        False, "--pep", help="Also export PEP metadata (config YAML + sample CSV)."
    ),
    resolver: str = RESOLVER,
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Resolve and report, but download nothing."
    ),
    quiet: bool = typer.Option(
        False, "--quiet", "-q", help="Suppress transfer progress output."
    ),
    layout: str = LAYOUT,
    verbose: bool = VERBOSE,
) -> None:
    """Download data files and metadata for each accession."""
    if convert is not None and convert != "fastq":
        error_console.print(
            f"[red]Unknown --convert target {convert!r}. Currently only "
            f"'fastq' is supported.[/red]"
        )
        raise typer.Exit(code=2)
    if convert == "fastq":
        # Resolve once, up front, so a missing toolkit fails before any
        # download starts rather than after hours of transfer.
        fasterq_dump, convert_env = discover(fasterq_dump, convert_env)
        if not fasterq_dump_available(fasterq_dump, convert_env):
            error_console.print(
                "[red]--convert fastq needs fasterq-dump, which could not be "
                "located.[/red]\n"
                "  Build it once with:  omixfetch setup\n"
                f"  Or point at a binary: export {PATH_ENV_VAR}=/path/to/fasterq-dump"
            )
            raise typer.Exit(code=2)
        where = f"environment {convert_env}" if convert_env else fasterq_dump
        console.print(f"[dim]fasterq-dump: {where}[/dim]")

    if layout not in LAYOUTS:
        error_console.print(
            f"[red]Unknown layout {layout!r}. Use one of: {', '.join(LAYOUTS)}.[/red]"
        )
        raise typer.Exit(code=2)

    limit = parse_size(max_size) if max_size else None
    want_reads = reads or reads_only
    datasets = _resolve_all(
        accessions,
        resolver=resolver,
        verbose=verbose,
        include_samples=samples,
        include_matrix=matrix,
        include_reads=want_reads,
        include_archives=include_archives,
        prefer=prefer,
        source=source,
        sra_flavour=sra_flavour,
    )
    if reads_only:
        kind = list(kind) + ["raw"] if kind else ["raw"]

    engine: object
    if backend == "auto":
        # s5cmd talks the S3 API directly and is usually faster against NCBI's
        # buckets, but only serves S3 URLs, so it is chosen only when every
        # selected file is there.
        on_s3 = bool(datasets) and all(
            s3_uri(f.url) for d in datasets for f in d.files
        )
        if on_s3 and s5cmd_available():
            engine_name = "s5cmd"
        elif aria2_available():
            engine_name = "aria2"
        else:
            engine_name = "python"
    else:
        engine_name = backend
    if engine_name == "s5cmd":
        engine = S5cmdBackend(concurrent=jobs, connections=connections, quiet=quiet)
    elif engine_name == "aspera":
        engine = AsperaBackend(quiet=quiet)
    elif engine_name == "aria2":
        engine = Aria2Backend(connections=connections, concurrent=jobs, quiet=quiet)
    elif engine_name == "python":
        if backend == "auto":
            error_console.print(
                "[yellow]aria2c not found; falling back to the built-in downloader. "
                "Install aria2 for multi-connection transfers.[/yellow]"
            )
        engine = RequestsBackend(concurrent=jobs, quiet=quiet)
    else:
        error_console.print(
            f"[red]Unknown backend {backend!r}. Use auto, aria2, s5cmd, "
            f"aspera or python.[/red]"
        )
        raise typer.Exit(code=2)

    for dataset in datasets:
        destination = outdir / dataset.accession
        if want_reads and not dataset.runs:
            error_console.print(
                f"[yellow]{dataset.accession}: ENA serves no public runs for this "
                f"project. Human data under controlled access (dbGaP) is the usual "
                f"reason; such reads need authorised access and cannot be fetched "
                f"here.[/yellow]"
            )
        selected = apply_filters(
            dataset.files, include=include, exclude=exclude, kinds=kind, max_size=limit
        )
        total = human_bytes(total_size(selected))
        console.print(
            f"[bold]{dataset.accession}[/bold]  {len(selected)} file(s), {total}"
            f"  ->  {destination}"
        )

        if metadata:
            write_metadata(dataset, destination / "metadata")
            if pep:
                write_pep(dataset, destination / "metadata" / "pep")

        apply_layout(selected, layout)

        if dry_run:
            for f in selected:
                console.print(f"    would fetch  {f.target}  ({human_bytes(f.size)})")
            continue

        if not selected:
            console.print("    nothing to download after filtering.")
            continue

        try:
            engine.download(selected, destination / "data")
        except OmixfetchError as exc:
            # Some environments have an aria2c that cannot complete a TLS
            # handshake (intercepting proxies, unusual TLS backends). When the
            # engine was chosen automatically, fall back rather than abort.
            if backend == "auto" and engine_name == "aria2":
                error_console.print(
                    f"[yellow]aria2c failed for {dataset.accession} ({exc}). "
                    f"Retrying with the built-in downloader.[/yellow]"
                )
                engine = RequestsBackend(concurrent=jobs, quiet=quiet)
                engine_name = "python"
                try:
                    engine.download(selected, destination / "data")
                    # aria2 keeps its input file for resume; the fallback
                    # succeeded, so there is nothing left to resume.
                    (destination / "data" / ".omixfetch.aria2.input").unlink(
                        missing_ok=True
                    )
                except OmixfetchError as fallback_exc:
                    error_console.print(
                        f"[red]error[/red] {dataset.accession}: {fallback_exc}"
                    )
                    raise typer.Exit(code=1) from fallback_exc
            else:
                error_console.print(f"[red]error[/red] {dataset.accession}: {exc}")
                raise typer.Exit(code=1) from exc
        if convert == "fastq":
            _convert_downloaded(
                destination / "data", threads, keep_sra, fasterq_dump, convert_env
            )

        console.print(f"    [green]complete[/green] -> {destination / 'data'}")


@app.command("fastq")
def fastq(
    accession: str = typer.Argument(..., help="A GEO series, e.g. GSE184708."),
    outdir: Path = typer.Option(Path("."), "--outdir", "-o", help="Destination directory."),
    route_order: str = typer.Option(
        "ena-fastq,sra,original",
        "--routes",
        help="Routes to try per run, in order. ena-fastq: ENA's own FASTQ; "
        "sra: NCBI .sra + fasterq-dump; original: submitted BAM + bamtofastq.",
    ),
    parallel_runs: int = typer.Option(
        2, "--parallel-runs", "-p", help="Runs processed at once."
    ),
    jobs: int = typer.Option(4, "--jobs", "-j", help="Files downloaded at once per run."),
    connections: int = typer.Option(8, "--connections", "-x", help="Connections per file."),
    threads: int = THREADS,
    include: list[str] = typer.Option(
        [], "--include", "-i", help="Only runs matching this glob, e.g. 'SRR160361*'."
    ),
    keep_intermediate: bool = typer.Option(
        False, "--keep-intermediate", help="Keep .sra/BAM files after conversion."
    ),
    disk_check: bool = typer.Option(
        True, "--disk-check/--no-disk-check", help="Refuse a run if free space looks short."
    ),
    fasterq_dump: str | None = FASTERQ_DUMP_PATH,
    convert_env: str | None = CONVERT_ENV,
    install: bool = typer.Option(
        True, "--install/--no-install",
        help="Install missing converters (SRA Toolkit, bamtofastq) automatically.",
    ),
    allow_incomplete: bool = typer.Option(
        False, "--allow-incomplete-sra",
        help="Use the .sra route even where it loses a read (10x R1). Rarely wanted.",
    ),
    dry_run: bool = typer.Option(False, "--dry-run", help="Plan only; download nothing."),
) -> None:
    """Fetch every run of a GEO series as verified FASTQ, autonomously.

    For each run the first working route is used, falling through on failure:
    ENA's own FASTQ, then NCBI's .sra converted with fasterq-dump, then the
    original submission converted with bamtofastq. Spot counts are checked
    against ENA's metadata. Progress is recorded, so rerunning the same
    command after an interruption resumes and skips finished runs.
    """
    import fnmatch
    import subprocess as _subprocess
    from concurrent.futures import ThreadPoolExecutor

    from .bundled import find_tool, install_bamtofastq, install_sra_toolkit
    from .convert import convert_sra as _convert_sra
    from .pipeline import (
        DONE,
        ROUTES,
        StateStore,
        Tools,
        is_complete,
        plan_runs,
        process_run,
        write_manifest,
    )
    from .probe import bam_header, is_tenx_bam

    order = tuple(r.strip() for r in route_order.split(",") if r.strip())
    unknown = [r for r in order if r not in ROUTES]
    if unknown:
        error_console.print(f"[red]Unknown route(s): {', '.join(unknown)}. "
                            f"Use: {', '.join(ROUTES)}.[/red]")
        raise typer.Exit(code=2)

    session = make_session()
    with console.status(f"Planning {accession} ..."):
        try:
            _dataset, plans = plan_runs(accession, session)
        except OmixfetchError as exc:
            error_console.print(f"[red]error[/red] {accession}: {exc}")
            raise typer.Exit(code=1) from exc

    if include:
        plans = [p for p in plans if any(fnmatch.fnmatch(p.run, g) for g in include)]
    if not plans:
        error_console.print(
            f"[yellow]{accession}: no public runs to fetch. Controlled-access "
            f"(dbGaP) data is the usual reason.[/yellow]"
        )
        raise typer.Exit(code=1)

    base = outdir / accession
    store = StateStore(base / ".omixfetch" / "state.json")

    # Converters: resolved once, up front, and installed if the plan needs them.
    binary, env_name = discover(fasterq_dump, convert_env)
    have_sra = fasterq_dump_available(binary, env_name)
    bamtofastq = find_tool("bamtofastq")

    def needs(route: str) -> bool:
        return any(
            route in p.available(order)
            and not (route == "sra" and p.incomplete_sra() and not allow_incomplete)
            for p in plans
        )

    if install and not dry_run:
        if needs("sra") and not have_sra:
            with console.status("Installing NCBI SRA Toolkit ..."):
                try:
                    binary, env_name = str(install_sra_toolkit(session)), None
                    have_sra = fasterq_dump_available(binary, None)
                except (OmixfetchError, OSError) as exc:
                    error_console.print(f"[yellow]SRA Toolkit install failed: {exc}[/yellow]")
        if needs("original") and not bamtofastq:
            with console.status("Installing 10x bamtofastq ..."):
                try:
                    bamtofastq = str(install_bamtofastq(session))
                except (OmixfetchError, OSError) as exc:
                    error_console.print(f"[yellow]bamtofastq install failed: {exc}[/yellow]")
    have_bam = bamtofastq is not None

    table = Table(title=f"{accession}: {len(plans)} run(s)")
    for col in ("Run", "Sample", "Spots", "Routes available", "Plan", "Status"):
        table.add_column(col)
    for plan in plans:
        usable = [r for r in plan.available(order)
                  if (r != "sra" or have_sra) and (r != "original" or have_bam)
                  and not (r == "sra" and plan.incomplete_sra() and not allow_incomplete)]
        sizes = ", ".join(
            f"{r} {human_bytes(sum(f.size or 0 for f in plan.routes[r]))}"
            for r in plan.available(order)
        )
        state = store.get(plan.run)
        table.add_row(
            plan.run, plan.sample or "",
            f"{plan.expected_spots:,}" if plan.expected_spots else "?",
            sizes or "none", usable[0] if usable else "[red]none usable[/red]",
            "done" if is_complete(state) else state.status,
        )
    console.print(table)
    console.print(
        f"[dim]fasterq-dump: {'environment ' + env_name if env_name else binary or 'not found'}"
        f" | bamtofastq: {bamtofastq or 'not found'}[/dim]"
    )
    flagged = [p for p in plans if p.incomplete_sra()]
    if flagged and not allow_incomplete:
        console.print(
            f"[yellow]{len(flagged)} run(s) skip the .sra route: "
            f"{flagged[0].incomplete_sra()}. Using the original submission.[/yellow]"
        )
    if dry_run:
        return

    log_path = base / ".omixfetch" / "pipeline.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_lock = __import__("threading").Lock()

    def log(message: str) -> None:
        console.print(message)
        with log_lock, open(log_path, "a") as handle:
            handle.write(time.strftime("%Y-%m-%d %H:%M:%S ") + message + "\n")

    def download(files, dest):
        # Prefer s5cmd for S3 objects, aria2 otherwise, Python as last resort.
        if all(s3_uri(f.url) for f in files) and s5cmd_available():
            S5cmdBackend(
                concurrent=jobs, connections=connections, quiet=True
            ).download(files, dest)
            return
        if aria2_available():
            try:
                Aria2Backend(
                    connections=connections, concurrent=jobs, quiet=True
                ).download(files, dest)
                return
            except OmixfetchError as exc:
                log(f"  aria2c failed ({str(exc)[:80]}); retrying with built-in downloader")
        RequestsBackend(concurrent=jobs, quiet=True).download(files, dest)

    def convert_bam(bam: Path, out: Path) -> list[Path]:
        target = out / "bamtofastq_out"
        proc = _subprocess.run(
            [bamtofastq, f"--nthreads={threads}", str(bam), str(target)],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            raise BackendError(f"bamtofastq failed: {(proc.stderr or '')[-300:]}")
        produced = []
        for path in sorted(target.rglob("*.fastq.gz")):
            # bamtofastq repeats filenames across flowcell folders; make unique.
            unique = path.with_name(f"{bam.stem}_{path.parent.name}_{path.name}")
            path.rename(unique)
            produced.append(unique)
        return produced

    def check_bam(url: str) -> None:
        # A few MB, read before committing to a ~100 GB download.
        try:
            header = bam_header(session, url)
        except (ValueError, OSError) as exc:
            raise BackendError(f"could not read the BAM header: {exc}") from exc
        if not is_tenx_bam(header):
            raise BackendError(
                "BAM has no 10x_bam_to_fastq header lines, so bamtofastq cannot "
                "rebuild reads from it (not a Cell Ranger BAM)"
            )

    tools = Tools(
        download=download,
        check_bam=check_bam,
        convert_sra=(
            (lambda sra, out: _convert_sra(sra, out, threads=threads,
                                           binary=binary, env_name=env_name))
            if have_sra else None
        ),
        convert_bam=convert_bam if have_bam else None,
        log=log,
    )

    log(f"{accession}: starting {len(plans)} run(s), {parallel_runs} at a time")
    with ThreadPoolExecutor(max_workers=max(1, parallel_runs)) as pool:
        list(pool.map(
            lambda plan: process_run(plan, base, store, tools, order=order,
                                     keep_intermediate=keep_intermediate,
                                     check_space=disk_check,
                                     allow_incomplete=allow_incomplete),
            plans,
        ))

    manifest = write_manifest(store, base / "fastq_manifest.tsv")
    states = [store.get(p.run) for p in plans]
    done = sum(1 for s in states if s.status == DONE)
    failed = [s for s in states if s.status != DONE]
    log(f"{accession}: {done}/{len(plans)} run(s) complete -> {base / 'fastq'}")
    console.print(f"manifest: {manifest}")
    for s in failed:
        error_console.print(f"  [red]{s.run}[/red] {s.status}: {(s.error or '')[:160]}")
    if failed:
        console.print("Rerun the same command to retry; finished runs are skipped.")
        raise typer.Exit(code=1)


@app.command("setup")
def setup(
    conda: bool = typer.Option(
        False, "--conda",
        help="Build a conda environment for sra-tools instead of using NCBI's "
        "prebuilt binaries.",
    ),
    env_name: str = typer.Option(
        DEFAULT_CONVERT_ENV, "--env-name", help="Name for the conda environment."
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Show what would be run, change nothing."
    ),
) -> None:
    """Install the converters: NCBI SRA Toolkit and 10x bamtofastq.

    By default both are fetched as prebuilt binaries into ~/.omixfetch/tools
    (override with OMIXFETCH_HOME), so no conda is involved -- which also
    avoids the osx-arm64 conflict between sra-tools and aria2. With --conda,
    sra-tools goes into its own conda environment instead.
    """
    from .bundled import (
        bamtofastq_asset,
        bin_dir,
        find_tool,
        install_bamtofastq,
        install_sra_toolkit,
        sra_toolkit_build,
    )

    if not conda:
        session = make_session()
        binary, _env = discover(None, None)
        if binary:
            console.print(f"fasterq-dump: {binary} (version {fasterq_dump_version(binary)})")
        else:
            console.print(f"Installing SRA Toolkit ({sra_toolkit_build()}) into {bin_dir()}")
            if not dry_run:
                path = install_sra_toolkit(session)
                console.print(f"  fasterq-dump {fasterq_dump_version(str(path))}: {path}")
        existing = find_tool("bamtofastq")
        if existing:
            console.print(f"bamtofastq: {existing}")
        else:
            console.print(f"Installing 10x bamtofastq ({bamtofastq_asset()}) into {bin_dir()}")
            if not dry_run:
                console.print(f"  bamtofastq: {install_bamtofastq(session)}")
        if not dry_run:
            console.print("[green]Done.[/green] fastq and --convert fastq will find these.")
        return

    _setup_conda(env_name, dry_run)


def _setup_conda(env_name: str, dry_run: bool) -> None:
    """Build the conversion environment that --convert fastq uses.

    sra-tools and aria2 cannot be solved into one conda environment on
    osx-arm64, so the toolkit lives in its own. This builds it once; after
    that --convert fastq finds it automatically, with no flags or exports.
    On Linux, where they coexist, an sra-tools already on PATH is used and
    this step is unnecessary.
    """
    binary, existing = discover(None, None)
    if binary:
        console.print(f"fasterq-dump already available: {binary}")
        console.print(f"  version {fasterq_dump_version(binary) or 'unknown'}")
        console.print("Nothing to do.")
        return
    if existing == env_name:
        console.print(f"Environment {env_name!r} already exists; verifying ...")
        ok = fasterq_dump_available(None, env_name)
        console.print("  [green]fasterq-dump works[/green]" if ok
                      else "  [red]fasterq-dump not runnable there[/red]")
        raise typer.Exit(code=0 if ok else 1)

    command = create_convert_env(env_name)
    console.print("Building the conversion environment:")
    console.print("  " + " ".join(command))
    if dry_run:
        return

    proc = subprocess.run(command)
    if proc.returncode != 0:
        error_console.print(
            f"[red]Environment creation failed (exit {proc.returncode}).[/red] "
            "NCBI also ships dependency-free binaries: "
            "https://github.com/ncbi/sra-tools/releases"
        )
        raise typer.Exit(code=1)

    if fasterq_dump_available(None, env_name):
        console.print(f"[green]Done.[/green] --convert fastq will use {env_name!r}.")
    else:
        error_console.print(
            f"[red]{env_name!r} was created but fasterq-dump does not run "
            "inside it.[/red]"
        )
        raise typer.Exit(code=1)


@app.command("benchmark")
def benchmark(
    accession: str = typer.Argument(..., help="A GEO accession to benchmark against."),
    sample: str = typer.Option(
        "100MB", "--sample", help="Bytes to sample from each source."
    ),
    resolver: str = RESOLVER,
) -> None:
    """Measure achievable throughput from each available source.

    Which mirror is fastest depends on the network path of this machine, so
    published figures are not transferable -- this samples a fixed number of
    bytes from one representative file per source and reports the rate, plus
    how long the whole dataset would take at that rate.
    """
    from .sources.ena import fetch_runs, project_accessions
    from .sources.ncbi import fetch_sdl, objects_to_files

    sample_bytes = parse_size(sample)
    session = make_session()

    with console.status(f"Resolving {accession} ..."):
        dataset = resolve(
            accession, session=session, resolver=resolver, include_reads=False
        )
        projects = project_accessions(dataset.metadata.relations)
        runs = fetch_runs(session, projects[0]) if projects else []

    if not runs:
        error_console.print(
            f"[red]{accession} has no public ENA runs to benchmark.[/red]"
        )
        raise typer.Exit(code=1)

    candidates: list[tuple[str, object]] = []
    ena_files = [f for run in runs for f in run.files]
    if ena_files:
        biggest = max(ena_files, key=lambda f: f.size or 0)
        candidates.append(("ENA HTTPS", biggest))

    objects = fetch_sdl(session, [r.accession for r in runs])
    totals = {}
    for flavour, label in (("sra", "NCBI S3 (.sra)"), ("original", "NCBI S3 (original)")):
        files = objects_to_files(objects, flavour=flavour)
        if files:
            candidates.append((label, max(files, key=lambda f: f.size or 0)))
            totals[label] = sum(f.size or 0 for f in files)
    totals["ENA HTTPS"] = sum(f.size or 0 for f in ena_files)

    console.print(f"[bold]{accession}[/bold]  {len(runs)} run(s)")
    for label, total in totals.items():
        console.print(f"  {label:<22} {human_bytes(total)} total")
    console.print(f"\nSampling {human_bytes(sample_bytes)} per source ...\n")

    results = []
    for label, f in candidates:
        with console.status(f"  {label} ..."):
            results.append(benchmark_source(label, f, session, sample_bytes))

    for line in format_table(results, None):
        console.print("  " + line)

    console.print()
    for r in results:
        total = totals.get(r.source)
        hours = r.projected_hours(total) if total else None
        if hours is not None:
            console.print(
                f"  {r.source:<22} full set ({human_bytes(total)}): {hours:.1f} h"
            )

    if not ascp_available():
        console.print(
            "\n  [dim]ascp not installed; ENA Aspera not benchmarked. "
            "conda install -c hcc aspera-cli[/dim]"
        )


def _convert_downloaded(
    data_dir: Path,
    threads: int,
    keep_sra: bool,
    binary: str | None = None,
    env_name: str | None = None,
) -> None:
    """Convert every .sra archive under *data_dir* to FASTQ, in place.

    Each archive is converted into its own directory so paired files stay
    beside the sample they came from.
    """
    archives = sorted(data_dir.rglob("*.sra"))
    if not archives:
        console.print("    no .sra archives to convert.")
        return

    for archive in archives:
        console.print(f"    converting  {archive.name} ...")
        try:
            produced = sra_to_fastq(
                archive,
                archive.parent,
                threads=threads,
                gzip_output=True,
                keep_sra=keep_sra,
                binary=binary,
                env_name=env_name,
            )
        except OmixfetchError as exc:
            error_console.print(f"[red]    conversion failed[/red]: {exc}")
            raise typer.Exit(code=1) from exc
        for path in produced:
            size = path.stat().st_size
            console.print(f"      wrote {path.name}  ({human_bytes(size)})")


if __name__ == "__main__":
    app()
