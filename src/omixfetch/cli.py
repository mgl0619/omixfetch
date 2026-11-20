"""Command-line interface for omixfetch."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

from . import __version__
from .backends import Aria2Backend, RequestsBackend, aria2_available, aria2_version
from .errors import OmixfetchError
from .filters import apply_filters, parse_size
from .http import make_session
from .manifest import human_bytes, total_size, write_metadata
from .models import Dataset
from .resolve import resolve

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Download public genomics datasets from NCBI GEO in parallel.",
)
console = Console()
error_console = Console(stderr=True)

ACCESSIONS = typer.Argument(..., help="GEO accessions (GSE, GSM or GPL).")


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
        engine = aria2_version() or "not installed"
        console.print(f"aria2c: {engine}")
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
def info(accessions: list[str] = ACCESSIONS) -> None:
    """Print study-level metadata for each accession."""
    for dataset in _resolve_all(accessions, include_samples=True, fetch_sizes=False):
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
    include: list[str] = typer.Option([], "--include", "-i", help="Keep files matching this glob."),
    exclude: list[str] = typer.Option([], "--exclude", "-e", help="Drop files matching this glob."),
    kind: list[str] = typer.Option([], "--kind", help="Keep only these file kinds."),
    max_size: Optional[str] = typer.Option(None, "--max-size", help="Drop files larger than this, e.g. 500MB."),
    samples: bool = typer.Option(True, "--samples/--no-samples", help="Include per-sample files."),
    sizes: bool = typer.Option(True, "--sizes/--no-sizes", help="Query file sizes (GEO needs one HEAD per file)."),
) -> None:
    """List the files an accession resolves to, without downloading."""
    limit = parse_size(max_size) if max_size else None
    for dataset in _resolve_all(accessions, include_samples=samples, fetch_sizes=sizes):
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
    samples: bool = typer.Option(True, "--samples/--no-samples", help="Include per-sample annotation."),
) -> None:
    """Fetch and write metadata only, downloading no data files."""
    for dataset in _resolve_all(accessions, include_samples=samples, fetch_sizes=False):
        directory = outdir / dataset.accession / "metadata"
        written = write_metadata(dataset, directory)
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
    outdir: Path = typer.Option(Path("."), "--outdir", help="Directory the aria2 manifest should target."),
    include: list[str] = typer.Option([], "--include", "-i"),
    exclude: list[str] = typer.Option([], "--exclude", "-e"),
    max_size: Optional[str] = typer.Option(None, "--max-size"),
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
    for dataset in _resolve_all(accessions, include_samples=True):
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
    connections: int = typer.Option(8, "--connections", "-x", help="Connections per file."),
    backend: str = typer.Option("auto", "--backend", "-b", help="auto, aria2 or python."),
    include: list[str] = typer.Option([], "--include", "-i", help="Keep files matching this glob."),
    exclude: list[str] = typer.Option([], "--exclude", "-e", help="Drop files matching this glob."),
    kind: list[str] = typer.Option([], "--kind", help="Keep only these file kinds."),
    max_size: Optional[str] = typer.Option(None, "--max-size", help="Skip files larger than this."),
    samples: bool = typer.Option(True, "--samples/--no-samples", help="Include per-sample files."),
    metadata: bool = typer.Option(True, "--metadata/--no-metadata", help="Write metadata alongside the data."),
    matrix: bool = typer.Option(False, "--matrix", help="Also fetch the GEO series matrix file."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Resolve and report, but download nothing."),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress transfer progress output."),
) -> None:
    """Download data files and metadata for each accession."""
    limit = parse_size(max_size) if max_size else None
    datasets = _resolve_all(accessions, include_samples=samples, include_matrix=matrix)

    engine: object
    if backend == "auto":
        engine_name = "aria2" if aria2_available() else "python"
    else:
        engine_name = backend
    if engine_name == "aria2":
        engine = Aria2Backend(connections=connections, concurrent=jobs, quiet=quiet)
    elif engine_name == "python":
        if backend == "auto":
            error_console.print(
                "[yellow]aria2c not found; falling back to the built-in downloader. "
                "Install aria2 for multi-connection transfers.[/yellow]"
            )
        engine = RequestsBackend(concurrent=jobs, quiet=quiet)
    else:
        error_console.print(f"[red]Unknown backend {backend!r}. Use auto, aria2 or python.[/red]")
        raise typer.Exit(code=2)

    for dataset in datasets:
        destination = outdir / dataset.accession
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

        if dry_run:
            for f in selected:
                console.print(f"    would fetch  {f.filename}  ({human_bytes(f.size)})")
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
                except OmixfetchError as fallback_exc:
                    error_console.print(
                        f"[red]error[/red] {dataset.accession}: {fallback_exc}"
                    )
                    raise typer.Exit(code=1) from fallback_exc
            else:
                error_console.print(f"[red]error[/red] {dataset.accession}: {exc}")
                raise typer.Exit(code=1) from exc
        console.print(f"    [green]complete[/green] -> {destination / 'data'}")


if __name__ == "__main__":
    app()
