"""s5cmd transfer backend for NCBI's public S3 mirrors.

NCBI serves both the original submitted files and the normalised ``.sra``
archives from public S3 buckets. s5cmd talks the S3 API directly and
parallelises across objects and multipart ranges, which is often faster than
generic HTTP clients against S3 -- particularly from a well-connected node.

Only NCBI's public buckets are supported, and they are requester-anonymous,
so transfers run with ``--no-sign-request`` and need no AWS credentials.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

from ..errors import BackendError
from ..models import RemoteFile

S5CMD = "s5cmd"

#: https://<bucket>.s3[.<region>].amazonaws.com/<key>
VIRTUAL_HOSTED_RE = re.compile(
    r"^https?://(?P<bucket>[^.]+)\.s3[.-][^/]*amazonaws\.com/(?P<key>.+)$",
    re.IGNORECASE,
)

#: https://s3[.<region>].amazonaws.com/<bucket>/<key>
PATH_STYLE_RE = re.compile(
    r"^https?://s3[.-][^/]*amazonaws\.com/(?P<bucket>[^/]+)/(?P<key>.+)$",
    re.IGNORECASE,
)


def s5cmd_available() -> bool:
    """Return True if an ``s5cmd`` executable is on PATH."""
    return shutil.which(S5CMD) is not None


def s5cmd_version() -> str | None:
    """Return the s5cmd version string, or None if unavailable."""
    if not s5cmd_available():
        return None
    try:
        proc = subprocess.run(
            [S5CMD, "version"], capture_output=True, text=True, timeout=30
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return (proc.stdout or proc.stderr or "").strip().splitlines()[0] if (
        proc.stdout or proc.stderr
    ) else None


def s3_uri(url: str) -> str | None:
    """Translate an HTTPS S3 URL into an ``s3://bucket/key`` URI.

    Returns None for anything that is not an Amazon S3 endpoint, so callers
    can fall back to an HTTP backend rather than failing.
    """
    for pattern in (VIRTUAL_HOSTED_RE, PATH_STYLE_RE):
        match = pattern.match(url.strip())
        if match:
            return f"s3://{match.group('bucket')}/{match.group('key')}"
    return None


def render_run_file(
    files: Sequence[RemoteFile], dest_dir: Path, concurrency: int = 16
) -> str:
    """Render an s5cmd command file.

    s5cmd runs a batch of commands from a file in parallel, which avoids one
    process per object and lets it schedule across the whole set.
    """
    lines = []
    for f in files:
        uri = s3_uri(f.url)
        if uri is None:
            continue
        target = dest_dir / Path(f.target)
        lines.append(f"cp --concurrency {concurrency} '{uri}' '{target}'")
    return "\n".join(lines) + ("\n" if lines else "")


class S5cmdBackend:
    """Download NCBI S3 objects with ``s5cmd``."""

    name = "s5cmd"

    def __init__(
        self,
        *,
        concurrent: int = 8,
        connections: int = 16,
        quiet: bool = False,
        **_: object,
    ) -> None:
        """
        Args:
            concurrent: objects transferred at once (``--numworkers``).
            connections: multipart parts per object (``cp --concurrency``).
        """
        if not s5cmd_available():
            raise BackendError(
                "s5cmd was not found on PATH. Install it with "
                "'conda install -c conda-forge s5cmd', or use --backend aria2."
            )
        self.concurrent = max(1, concurrent)
        self.connections = max(1, connections)
        self.quiet = quiet

    def build_command(self, run_file: Path) -> list[str]:
        """Return the s5cmd argument vector for a batch command file."""
        return [
            S5CMD,
            "--no-sign-request",
            "--numworkers",
            str(self.concurrent),
            *(["--log", "error"] if self.quiet else []),
            "run",
            str(run_file),
        ]

    def download(self, files: Sequence[RemoteFile], dest_dir: Path) -> Path:
        """Download *files* into *dest_dir*.

        Raises:
            BackendError: if any file is not an S3 URL, or s5cmd fails.
        """
        dest_dir = Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)
        if not files:
            return dest_dir

        non_s3 = [f.filename for f in files if s3_uri(f.url) is None]
        if non_s3:
            raise BackendError(
                f"{len(non_s3)} file(s) are not on S3 (first: {non_s3[0]}). "
                f"s5cmd only serves NCBI's S3 mirrors -- use --source ncbi, "
                f"or --backend aria2 for these."
            )

        for f in files:
            (dest_dir / Path(f.target)).parent.mkdir(parents=True, exist_ok=True)

        run_file = dest_dir / ".omixfetch.s5cmd.txt"
        run_file.write_text(render_run_file(files, dest_dir, self.connections))

        proc = subprocess.run(self.build_command(run_file), check=False)
        if proc.returncode != 0:
            raise BackendError(
                f"s5cmd exited with status {proc.returncode}. The command file "
                f"was kept at {run_file}."
            )
        run_file.unlink(missing_ok=True)
        return dest_dir
