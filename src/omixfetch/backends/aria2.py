"""aria2c transfer backend.

aria2c is the default engine because it already solves the parts of parallel
downloading that are easy to get subtly wrong: multi-connection range
requests against a single large file, correct resumption of a partially
transferred multi-connection download, per-file checksum verification and
retry with backoff.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from collections.abc import Iterable, Sequence
from pathlib import Path

from ..errors import BackendError
from ..models import RemoteFile


def aria2_available() -> bool:
    """Return True if an ``aria2c`` executable is on PATH."""
    return shutil.which("aria2c") is not None


def aria2_version() -> str | None:
    """Return the aria2c version string, or None if unavailable."""
    exe = shutil.which("aria2c")
    if exe is None:
        return None
    try:
        proc = subprocess.run(
            [exe, "--version"], capture_output=True, text=True, timeout=30, check=False
        )
    except (OSError, subprocess.SubprocessError):
        return None
    first = proc.stdout.splitlines()[0] if proc.stdout else ""
    return first.strip() or None


def ca_bundle() -> str | None:
    """Return a CA bundle path for aria2c, or None to use its own default.

    aria2c builds against several TLS backends and does not always find the
    system trust store -- notably inside conda environments and containers,
    where the handshake fails with "Invalid certificate format". requests
    ships certifi, so a usable bundle is almost always on hand.
    """
    for key in ("SSL_CERT_FILE", "REQUESTS_CA_BUNDLE", "CURL_CA_BUNDLE"):
        value = os.environ.get(key)
        if value and os.path.exists(value):
            return value
    try:
        import certifi
    except ImportError:
        return None
    path = certifi.where()
    return path if os.path.exists(path) else None


def subprocess_env() -> dict[str, str]:
    """Return an environment safe for aria2c.

    aria2c understands HTTP proxies only. A SOCKS URL in ``all_proxy`` -- which
    many sandboxes and corporate networks set -- makes it abort at start-up
    with "unrecognized proxy format", so such values are dropped rather than
    passed through. HTTP proxy variables are left untouched.
    """
    env = dict(os.environ)
    for key in ("all_proxy", "ALL_PROXY"):
        value = env.get(key, "")
        if value and not value.lower().startswith(("http://", "https://")):
            env.pop(key, None)
    return env


def render_input_file(files: Sequence[RemoteFile], dest_dir: Path) -> str:
    """Render an aria2c input file for *files* landing under *dest_dir*.

    The input-file format puts per-download options on indented lines after
    the URL, which lets every file carry its own output name and, when the
    archive publishes one, its own checksum.
    """
    lines: list[str] = []
    for f in files:
        target = Path(f.target)
        lines.append(f.url)
        # aria2 creates the directory, so grouping is expressed through dir=.
        directory = dest_dir / target.parent if target.parent != Path(".") else dest_dir
        lines.append(f"  dir={directory}")
        lines.append(f"  out={target.name}")
        if f.md5:
            lines.append(f"  checksum=md5={f.md5}")
    return "\n".join(lines) + ("\n" if lines else "")


class Aria2Backend:
    """Download files by delegating to an ``aria2c`` subprocess."""

    name = "aria2"

    def __init__(
        self,
        *,
        connections: int = 8,
        concurrent: int = 4,
        max_tries: int = 5,
        retry_wait: int = 5,
        quiet: bool = False,
        extra_args: Iterable[str] = (),
    ) -> None:
        """
        Args:
            connections: range connections opened per file (aria2 ``-x``/``-s``).
            concurrent: files transferred at once (aria2 ``-j``).
            max_tries: attempts per file before giving up.
            retry_wait: seconds between attempts.
            quiet: suppress aria2's progress output.
            extra_args: additional raw aria2c arguments.
        """
        if not aria2_available():
            raise BackendError(
                "aria2c was not found on PATH. Install it (macOS: 'brew install aria2', "
                "Debian/Ubuntu: 'apt install aria2', conda: 'conda install -c conda-forge "
                "aria2'), or rerun with '--backend python' to use the slower built-in "
                "downloader."
            )
        self.connections = max(1, connections)
        self.concurrent = max(1, concurrent)
        self.max_tries = max_tries
        self.retry_wait = retry_wait
        self.quiet = quiet
        self.extra_args = list(extra_args)

    def build_command(self, input_path: Path) -> list[str]:
        """Return the full aria2c argument vector for *input_path*."""
        return [
            "aria2c",
            f"--input-file={input_path}",
            f"--max-concurrent-downloads={self.concurrent}",
            f"--max-connection-per-server={self.connections}",
            f"--split={self.connections}",
            "--min-split-size=10M",
            "--continue=true",
            "--auto-file-renaming=false",
            "--allow-overwrite=false",
            "--conditional-get=true",
            f"--max-tries={self.max_tries}",
            f"--retry-wait={self.retry_wait}",
            "--timeout=60",
            "--connect-timeout=30",
            "--file-allocation=none",
            *([f"--ca-certificate={ca}"] if (ca := ca_bundle()) else []),
            "--summary-interval=0",
            f"--console-log-level={'warn' if self.quiet else 'notice'}",
            *(["--quiet=true"] if self.quiet else []),
            *self.extra_args,
        ]

    def download(self, files: Sequence[RemoteFile], dest_dir: Path) -> Path:
        """Download *files* into *dest_dir*.

        Returns:
            The directory the files were written to.

        Raises:
            BackendError: if aria2c exits non-zero.
        """
        dest_dir = Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)
        if not files:
            return dest_dir

        input_path = dest_dir / ".omixfetch.aria2.input"
        input_path.write_text(render_input_file(files, dest_dir))

        command = self.build_command(input_path)
        proc = subprocess.run(command, check=False, env=subprocess_env())
        if proc.returncode != 0:
            raise BackendError(
                f"aria2c exited with status {proc.returncode}. The input file was kept "
                f"at {input_path} -- rerun the same command to resume."
            )
        input_path.unlink(missing_ok=True)
        return dest_dir
