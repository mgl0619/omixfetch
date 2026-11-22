"""Convert downloaded ``.sra`` archives to FASTQ.

Downloading and converting are two steps that almost always go together, so
omixfetch runs them as one. The conversion itself is delegated to
``fasterq-dump`` from the SRA Toolkit -- it is the reference implementation
and reimplementing it would be a mistake.

The toolkit is NOT a Python dependency and pip will not install it. Install
it separately:

    mamba install -c bioconda sra-tools
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

from .errors import BackendError

FASTERQ_DUMP = "fasterq-dump"


def fasterq_dump_available() -> bool:
    """Return True if ``fasterq-dump`` is on PATH."""
    return shutil.which(FASTERQ_DUMP) is not None


def fasterq_dump_version() -> str | None:
    """Return the installed fasterq-dump version, or None."""
    if not fasterq_dump_available():
        return None
    try:
        proc = subprocess.run(
            [FASTERQ_DUMP, "--version"], capture_output=True, text=True, timeout=60
        )
    except (OSError, subprocess.SubprocessError):
        return None
    for line in (proc.stdout or "").splitlines():
        if ":" in line:
            return line.split(":", 1)[1].strip()
    return None


def build_command(
    sra_path: Path,
    outdir: Path,
    *,
    threads: int = 4,
    split: str = "3",
    temp_dir: Path | None = None,
) -> list[str]:
    """Return the fasterq-dump argument vector.

    Args:
        split: ``"3"`` writes ``_1``/``_2`` for paired runs and a single file
            for unpaired ones, which is what almost everyone wants;
            ``"files"`` always splits; ``"spot"`` writes one interleaved file.

    fasterq-dump needs scratch space roughly the size of the output, so the
    temporary directory is pointed at the destination rather than ``/tmp``,
    which is often small on a cluster node.
    """
    flag = {"3": "--split-3", "files": "--split-files", "spot": "--split-spot"}
    if split not in flag:
        raise ValueError(f"Unknown split mode {split!r}. Use 3, files or spot.")

    return [
        FASTERQ_DUMP,
        str(sra_path),
        "--outdir",
        str(outdir),
        "--temp",
        str(temp_dir or outdir),
        "--threads",
        str(max(1, threads)),
        flag[split],
        "--skip-technical",
        "--force",
    ]


def compress(path: Path, remove_original: bool = True) -> Path:
    """gzip *path* in place, returning the compressed path.

    fasterq-dump writes plain FASTQ; almost every downstream tool wants it
    gzipped, and uncompressed FASTQ is roughly three times the size.
    """
    target = path.with_suffix(path.suffix + ".gz")
    with open(path, "rb") as src, gzip.open(target, "wb") as dst:
        shutil.copyfileobj(src, dst, length=1 << 20)
    if remove_original:
        path.unlink()
    return target


def sra_to_fastq(
    sra_path: Path,
    outdir: Path,
    *,
    threads: int = 4,
    split: str = "3",
    gzip_output: bool = True,
    keep_sra: bool = False,
) -> list[Path]:
    """Convert one ``.sra`` archive to FASTQ files.

    Returns:
        The FASTQ files written, gzipped unless disabled.

    Raises:
        BackendError: if the toolkit is missing or the conversion fails.
    """
    if not fasterq_dump_available():
        raise BackendError(
            "fasterq-dump was not found on PATH. omixfetch can download .sra "
            "archives without the SRA Toolkit, but converting them to FASTQ "
            "needs it. Install with: mamba install -c bioconda sra-tools"
        )

    sra_path = Path(sra_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    before = set(outdir.glob("*.fastq"))

    command = build_command(sra_path, outdir, threads=threads, split=split)
    proc = subprocess.run(command, capture_output=True, text=True)
    if proc.returncode != 0:
        raise BackendError(
            f"fasterq-dump failed for {sra_path.name} "
            f"(exit {proc.returncode}): {(proc.stderr or '').strip()[:400]}"
        )

    produced = sorted(set(outdir.glob("*.fastq")) - before)
    if not produced:
        raise BackendError(
            f"fasterq-dump reported success for {sra_path.name} but wrote no "
            f"FASTQ files into {outdir}."
        )

    if gzip_output:
        produced = [compress(p) for p in produced]
    if not keep_sra:
        sra_path.unlink(missing_ok=True)
    return produced
