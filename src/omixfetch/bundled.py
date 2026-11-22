"""Install the external converters omixfetch needs into its own directory.

Two tools turn downloaded archives into FASTQ:

* ``fasterq-dump`` from NCBI's SRA Toolkit, for ``.sra`` archives.
* 10x Genomics ``bamtofastq``, for Cell Ranger BAMs. This is the only route
  that recovers cell barcodes and UMIs for 10x submissions: the ``.sra`` of
  a 10x BAM keeps only the cDNA read (GSE184708 has 98 bases per spot and no
  mates), while the BAM keeps barcodes in its tags.

Both vendors publish self-contained prebuilt binaries. Fetching those avoids
conda entirely -- including the osx-arm64 conflict between sra-tools and
aria2, which cannot be solved into one conda environment.

Tools are installed under ``$OMIXFETCH_HOME/tools`` (default
``~/.omixfetch/tools``). On a cluster with a small home quota, point
OMIXFETCH_HOME somewhere larger; the toolkit is ~200 MB unpacked.
"""

from __future__ import annotations

import os
import platform
import shutil
import stat
import tarfile
from pathlib import Path

import requests

from .errors import BackendError

SRA_BASE = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current"
BAMTOFASTQ_RELEASE = "https://api.github.com/repos/10XGenomics/bamtofastq/releases/latest"

HOME_VAR = "OMIXFETCH_HOME"


def home() -> Path:
    """Root directory for omixfetch's own files."""
    return Path(os.environ.get(HOME_VAR) or Path.home() / ".omixfetch")


def tools_dir() -> Path:
    return home() / "tools"


def bin_dir() -> Path:
    """Directory holding links to every installed tool."""
    return tools_dir() / "bin"


def _linux_flavour() -> str:
    """NCBI ships separate builds for Debian-like and RHEL-like systems."""
    try:
        text = Path("/etc/os-release").read_text().lower()
    except OSError:
        return "ubuntu64"
    rhel = ("rhel", "centos", "rocky", "almalinux", "alma", "fedora", "ol")
    for line in text.splitlines():
        if line.startswith(("id=", "id_like=")):
            values = line.split("=", 1)[1].strip().strip('"').split()
            if any(v in rhel for v in values):
                return "alma_linux64"
    return "ubuntu64"


def sra_toolkit_build(system: str | None = None, machine: str | None = None) -> str:
    """Return NCBI's build name for this platform.

    Raises:
        BackendError: on a platform NCBI does not build for.
    """
    system = (system or platform.system()).lower()
    machine = (machine or platform.machine()).lower()
    if system == "darwin":
        return "mac-arm64" if machine in ("arm64", "aarch64") else "mac64"
    if system == "linux" and machine in ("x86_64", "amd64"):
        return _linux_flavour()
    raise BackendError(
        f"NCBI publishes no prebuilt SRA Toolkit for {system}/{machine}. "
        "Install it another way and point omixfetch at it with --fasterq-dump."
    )


def bamtofastq_asset(system: str | None = None) -> str:
    """Return the release asset name for this platform."""
    system = (system or platform.system()).lower()
    if system == "darwin":
        return "bamtofastq_macos"
    if system == "linux":
        return "bamtofastq_linux"
    raise BackendError(f"10x publishes no bamtofastq build for {system}.")


def _download(session: requests.Session, url: str, target: Path) -> Path:
    """Stream *url* to *target*, verifying the declared length."""
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    with session.get(url, stream=True, timeout=300) as response:
        response.raise_for_status()
        expected = response.headers.get("Content-Length")
        with open(partial, "wb") as handle:
            for chunk in response.iter_content(1 << 20):
                handle.write(chunk)
    if expected and partial.stat().st_size != int(expected):
        partial.unlink(missing_ok=True)
        raise BackendError(f"incomplete download of {url}")
    partial.replace(target)
    return target


def _link(binary: Path, name: str) -> Path:
    """Expose *binary* as ``bin_dir()/name``."""
    bin_dir().mkdir(parents=True, exist_ok=True)
    link = bin_dir() / name
    if link.is_symlink() or link.exists():
        link.unlink()
    # Absolute target, so the link survives a change of working directory.
    link.symlink_to(Path(binary).resolve())
    return link


def _safe_extract(archive: Path, dest: Path) -> None:
    """Extract a tarball, refusing members that escape *dest*."""
    dest = dest.resolve()
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            target = (dest / member.name).resolve()
            if dest not in target.parents and target != dest:
                raise BackendError(f"unsafe path in archive: {member.name}")
        tar.extractall(dest)


def install_sra_toolkit(session: requests.Session | None = None) -> Path:
    """Download and unpack NCBI's prebuilt SRA Toolkit; return fasterq-dump."""
    session = session or requests.Session()
    build = sra_toolkit_build()
    url = f"{SRA_BASE}/sratoolkit.current-{build}.tar.gz"
    root = tools_dir() / "sratoolkit"
    archive = _download(session, url, tools_dir() / "downloads" / f"sratoolkit-{build}.tar.gz")

    staging = root.with_name("sratoolkit.staging")
    shutil.rmtree(staging, ignore_errors=True)
    staging.mkdir(parents=True)
    _safe_extract(archive, staging)
    shutil.rmtree(root, ignore_errors=True)
    staging.replace(root)
    archive.unlink(missing_ok=True)

    binary = next(root.glob("*/bin/fasterq-dump"), None)
    if binary is None:
        raise BackendError("SRA Toolkit unpacked but fasterq-dump was not found in it.")
    for tool in ("fasterq-dump", "prefetch", "vdb-config", "sam-dump"):
        candidate = binary.parent / tool
        if candidate.exists():
            _link(candidate, tool)
    return bin_dir() / "fasterq-dump"


def install_bamtofastq(session: requests.Session | None = None) -> Path:
    """Download 10x's prebuilt bamtofastq; return its path."""
    session = session or requests.Session()
    asset = bamtofastq_asset()
    release = session.get(BAMTOFASTQ_RELEASE, timeout=60)
    release.raise_for_status()
    url = next(
        (a["browser_download_url"] for a in release.json().get("assets", [])
         if a.get("name") == asset),
        None,
    )
    if url is None:
        raise BackendError(f"the latest bamtofastq release has no {asset} asset")
    target = _download(session, url, tools_dir() / "bamtofastq" / asset)
    target.chmod(target.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return _link(target, "bamtofastq")


def bundled(name: str) -> str | None:
    """Return the bundled copy of *name* if it is installed and executable."""
    path = bin_dir() / name
    if path.exists() and os.access(path, os.X_OK):
        return str(path)
    return None


def find_tool(name: str) -> str | None:
    """Locate *name*: PATH first, then omixfetch's bundled copy."""
    return shutil.which(name) or bundled(name)
