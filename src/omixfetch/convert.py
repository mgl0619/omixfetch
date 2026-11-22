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
import os
import shutil
import subprocess
from pathlib import Path

from .errors import BackendError

FASTERQ_DUMP = "fasterq-dump"

#: Point this at a fasterq-dump binary, or at the directory holding it, to use
#: a toolkit installed outside omixfetch's own environment.
PATH_ENV_VAR = "OMIXFETCH_FASTERQ_DUMP"

#: Signatures of a conversion that failed because the archive is
#: reference-compressed and the reference could not be fetched.
_REFERENCE_ERRORS = ("as a RefSeq or as a WGS", "rcVFS", "cannot open reference")


def resolve_fasterq_dump(explicit: str | None = None) -> str | None:
    """Locate a usable ``fasterq-dump``, or return None.

    Resolution order: an explicitly supplied path, then
    ``OMIXFETCH_FASTERQ_DUMP``, then PATH.

    The toolkit deliberately does not have to live in omixfetch's own
    environment. On osx-arm64 sra-tools and aria2 cannot be solved together
    -- installing sra-tools removes aria2 -- so the converter commonly lives
    in a separate environment, a container, or a standalone NCBI tarball.
    """
    for candidate in (explicit, os.environ.get(PATH_ENV_VAR)):
        if not candidate:
            continue
        path = Path(candidate)
        if path.is_dir():
            path = path / FASTERQ_DUMP
        if path.exists() and os.access(path, os.X_OK):
            return str(path)
    return shutil.which(FASTERQ_DUMP)


def fasterq_dump_available(explicit: str | None = None) -> bool:
    """Return True if a fasterq-dump binary can be located."""
    return resolve_fasterq_dump(explicit) is not None


def fasterq_dump_version(explicit: str | None = None) -> str | None:
    """Return the fasterq-dump version, or None if it cannot be located."""
    binary = resolve_fasterq_dump(explicit)
    if binary is None:
        return None
    try:
        proc = subprocess.run(
            [binary, "--version"], capture_output=True, text=True, timeout=60
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
    binary: str | None = None,
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
        binary or FASTERQ_DUMP,
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
    binary: str | None = None,
) -> list[Path]:
    """Convert one ``.sra`` archive to FASTQ files.

    Returns:
        The FASTQ files written, gzipped unless disabled.

    Raises:
        BackendError: if the toolkit is missing or the conversion fails.
    """
    resolved = resolve_fasterq_dump(binary)
    if resolved is None:
        raise BackendError(
            "fasterq-dump could not be located. Downloading .sra archives does "
            "not need the SRA Toolkit, but converting them to FASTQ does.\n"
            "  - install it separately: mamba env create -f environment-convert.yml\n"
            f"  - or point omixfetch at an existing binary: "
            f"export {PATH_ENV_VAR}=/path/to/fasterq-dump\n"
            "  - or use --fasterq-dump /path/to/fasterq-dump\n"
            "It does not have to live in omixfetch's environment; on osx-arm64 "
            "sra-tools and aria2 cannot be installed together."
        )

    sra_path = Path(sra_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    before = set(outdir.glob("*.fastq"))

    command = build_command(
        sra_path, outdir, threads=threads, split=split, binary=resolved
    )
    proc = subprocess.run(command, capture_output=True, text=True)
    if proc.returncode != 0:
        stderr = (proc.stderr or "").strip()
        if any(sig in stderr for sig in _REFERENCE_ERRORS):
            raise BackendError(
                f"fasterq-dump could not reconstruct reads from {sra_path.name}: "
                f"this archive is reference-compressed, and the reference "
                f"sequences it needs could not be fetched from NCBI.\n"
                f"  - aligned submissions (BAM/CRAM) are stored this way, which "
                f"is part of why the .sra is smaller than the original\n"
                f"  - conversion therefore needs outbound access to NCBI's "
                f"reference service; check it is reachable and that "
                f"`vdb-config --report-cloud-identity` is sane\n"
                f"  - or avoid the reference entirely by fetching the "
                f"submitter's own file instead: --sra-flavour original\n"
                f"fasterq-dump said: {stderr[:300]}"
            )
        raise BackendError(
            f"fasterq-dump failed for {sra_path.name} "
            f"(exit {proc.returncode}): {stderr[:400]}"
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
