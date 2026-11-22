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
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .errors import BackendError

FASTERQ_DUMP = "fasterq-dump"

#: Point this at a fasterq-dump binary, or at the directory holding it, to use
#: a toolkit installed outside omixfetch's own environment.
PATH_ENV_VAR = "OMIXFETCH_FASTERQ_DUMP"

#: Environment variable naming a conda environment to run the toolkit in.
ENV_NAME_VAR = "OMIXFETCH_CONVERT_ENV"

#: Name of the conversion environment `omixfetch setup` creates, and the one
#: looked for automatically when fasterq-dump is not on PATH.
DEFAULT_CONVERT_ENV = "omixfetch-convert"

#: Packages installed into the conversion environment. sra-tools 3.4+ avoids
#: the perl-xml-libxml dependency that has no osx-arm64 build.
CONVERT_PACKAGES = ("sra-tools>=3.4",)
CONVERT_CHANNELS = ("conda-forge", "bioconda")

#: Launchers that can run a command inside another environment, preferred
#: order; micromamba and mamba start faster than conda.
_RUNNERS = ("micromamba", "mamba", "conda")

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


def find_runner() -> str | None:
    """Return a launcher able to run a command inside another environment."""
    for runner in _RUNNERS:
        path = shutil.which(runner)
        if path:
            return path
    return None


def env_prefix(env_name: str | None) -> list[str]:
    """Return the argv prefix that runs a command inside *env_name*.

    ``run -n <env>`` executes in another environment without activating it,
    which is how one omixfetch invocation can download in one environment and
    convert in another -- necessary because sra-tools and aria2 cannot share
    an environment on osx-arm64.

    Raises:
        BackendError: if no conda-family launcher is available.
    """
    if not env_name:
        return []
    runner = find_runner()
    if runner is None:
        raise BackendError(
            f"Cannot run the toolkit in environment {env_name!r}: none of "
            f"{', '.join(_RUNNERS)} is on PATH. Point at the binary directly "
            f"instead, with --fasterq-dump /path/to/fasterq-dump."
        )
    return [runner, "run", "-n", env_name]


def env_exists(env_name: str) -> bool:
    """Return True if a conda environment called *env_name* exists."""
    runner = find_runner()
    if runner is None:
        return False
    try:
        proc = subprocess.run(
            [runner, "env", "list"], capture_output=True, text=True, timeout=60
        )
    except (OSError, subprocess.SubprocessError):
        return False
    for line in (proc.stdout or "").splitlines():
        fields = line.split()
        if not fields or line.startswith("#"):
            continue
        if fields[0] == env_name or line.rstrip().endswith(f"/{env_name}"):
            return True
    return False


def discover(
    explicit: str | None = None, env_name: str | None = None
) -> tuple[str | None, str | None]:
    """Find fasterq-dump without the user having to configure anything.

    Returns ``(binary, env_name)``: exactly one is set when something usable
    was found, both are None otherwise. Resolution order:

    1. an explicit ``--fasterq-dump`` path or ``OMIXFETCH_FASTERQ_DUMP``
    2. an explicit ``--convert-env`` or ``OMIXFETCH_CONVERT_ENV``
    3. ``fasterq-dump`` on this PATH
    4. the copy ``omixfetch setup`` installs under ``~/.omixfetch/tools``
    5. an ``omixfetch-convert`` conda environment, if one was built
    """
    for candidate in (explicit, os.environ.get(PATH_ENV_VAR)):
        if candidate:
            binary = resolve_fasterq_dump(candidate)
            if binary:
                return binary, None

    named = env_name or os.environ.get(ENV_NAME_VAR)
    if named:
        return None, named

    on_path = shutil.which(FASTERQ_DUMP)
    if on_path:
        return on_path, None

    from .bundled import bundled

    installed = bundled(FASTERQ_DUMP)
    if installed:
        return installed, None

    if env_exists(DEFAULT_CONVERT_ENV):
        return None, DEFAULT_CONVERT_ENV
    return None, None


def create_convert_env(env_name: str = DEFAULT_CONVERT_ENV) -> list[str]:
    """Return the command that builds the conversion environment.

    Kept separate from execution so it can be shown in a dry run and tested.
    """
    runner = find_runner()
    if runner is None:
        raise BackendError(
            "Building the conversion environment needs micromamba, mamba or "
            "conda on PATH."
        )
    channels = [arg for ch in CONVERT_CHANNELS for arg in ("-c", ch)]
    return [runner, "create", "-y", "-n", env_name, *channels, *CONVERT_PACKAGES]


def fasterq_dump_available(
    explicit: str | None = None, env_name: str | None = None
) -> bool:
    """Return True if fasterq-dump can be located, or reached via *env_name*.

    When an environment is named the binary lives elsewhere by design, so it
    is confirmed by asking that environment rather than searching this PATH.
    """
    if env_name:
        try:
            prefix = env_prefix(env_name)
        except BackendError:
            return False
        try:
            proc = subprocess.run(
                [*prefix, FASTERQ_DUMP, "--version"],
                capture_output=True,
                text=True,
                timeout=180,
            )
        except (OSError, subprocess.SubprocessError):
            return False
        return proc.returncode == 0
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
    env_name: str | None = None,
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
        *env_prefix(env_name),
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


@dataclass
class ConversionResult:
    """Outcome of one fasterq-dump run."""

    files: list[Path]
    #: Spots fasterq-dump reports reading; compared against archive metadata.
    spots_read: int | None = None
    reads_written: int | None = None


_STAT_RE = re.compile(r"^\s*(spots read|reads written)\s*:\s*([\d,]+)", re.M)


def parse_stats(text: str) -> tuple[int | None, int | None]:
    """Pull 'spots read' and 'reads written' out of fasterq-dump's report."""
    found = {k: int(v.replace(",", "")) for k, v in _STAT_RE.findall(text or "")}
    return found.get("spots read"), found.get("reads written")


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
    env_name: str | None = None,
) -> list[Path]:
    """Convert one ``.sra`` archive to FASTQ files; see :func:`convert_sra`."""
    return convert_sra(
        sra_path,
        outdir,
        threads=threads,
        split=split,
        gzip_output=gzip_output,
        keep_sra=keep_sra,
        binary=binary,
        env_name=env_name,
    ).files


def convert_sra(
    sra_path: Path,
    outdir: Path,
    *,
    threads: int = 4,
    split: str = "3",
    gzip_output: bool = True,
    keep_sra: bool = False,
    binary: str | None = None,
    env_name: str | None = None,
) -> ConversionResult:
    """Convert one ``.sra`` archive to FASTQ files.

    Returns:
        The FASTQ files written, gzipped unless disabled.

    Raises:
        BackendError: if the toolkit is missing or the conversion fails.
    """
    # With an environment named, the binary is found there, not on this PATH.
    resolved = FASTERQ_DUMP if env_name else resolve_fasterq_dump(binary)
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
        sra_path,
        outdir,
        threads=threads,
        split=split,
        binary=resolved,
        env_name=env_name,
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

    spots, written = parse_stats((proc.stdout or "") + "\n" + (proc.stderr or ""))
    if gzip_output:
        produced = [compress(p) for p in produced]
    if not keep_sra:
        sra_path.unlink(missing_ok=True)
    return ConversionResult(files=produced, spots_read=spots, reads_written=written)
