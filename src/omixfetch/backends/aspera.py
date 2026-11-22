"""Aspera (``ascp``) transfer backend for ENA.

ENA publishes an Aspera endpoint alongside every FTP one. Aspera uses a UDP
protocol that is frequently much faster than FTP or HTTPS over long, lossy
paths -- which is exactly the situation when pulling terabytes from EBI to a
cluster elsewhere.

Two things have to be true for it to work, and neither can be assumed:

* ``ascp`` must be installed (``conda install -c hcc aspera-cli``).
* Outbound UDP on port 33001 must be permitted. Many institutional firewalls
  block it, in which case transfers fail rather than falling back, so the
  backend reports the problem clearly instead of hanging.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

from ..errors import BackendError
from ..models import RemoteFile

ASCP = "ascp"

#: ENA's public Aspera account.
ENA_USER = "era-fasp"

#: Key filenames shipped with aspera-cli / Aspera Connect.
KEY_NAMES = ("asperaweb_id_dsa.openssh", "aspera_bypass_dsa.pem")


def ascp_available() -> bool:
    """Return True if an ``ascp`` executable is on PATH."""
    return shutil.which(ASCP) is not None


def find_key() -> str | None:
    """Locate the Aspera private key that ships with the client.

    aspera-cli installs it under the conda prefix; Aspera Connect puts it in
    the user's home directory. ``ASPERA_SCP_KEY`` overrides both.
    """
    override = os.environ.get("ASPERA_SCP_KEY")
    if override and os.path.exists(override):
        return override

    roots = [
        os.environ.get("CONDA_PREFIX", ""),
        os.path.expanduser("~/.aspera/connect"),
        "/opt/aspera",
        "/usr/local/aspera",
    ]
    for root in filter(None, roots):
        for base, _dirs, names in os.walk(root):
            for name in KEY_NAMES:
                if name in names:
                    return os.path.join(base, name)
    return None


def remote_spec(file: RemoteFile) -> str:
    """Return the ``user@host:/path`` specification for *file*.

    ENA reports endpoints as ``fasp.sra.ebi.ac.uk:/vol1/...``; the account
    name has to be prepended.
    """
    endpoint = (file.aspera or "").strip()
    if not endpoint:
        raise BackendError(f"{file.filename}: no Aspera endpoint published.")
    if "@" in endpoint.split(":", 1)[0]:
        return endpoint
    return f"{ENA_USER}@{endpoint}"


class AsperaBackend:
    """Transfer files with ``ascp``."""

    name = "aspera"

    def __init__(
        self,
        *,
        target_rate: str = "300m",
        key: str | None = None,
        port: int = 33001,
        quiet: bool = False,
        **_: object,
    ) -> None:
        """
        Args:
            target_rate: requested rate, e.g. ``"300m"``. Aspera does not
                probe for available bandwidth, so an unrealistically high
                figure simply wastes effort while too low a one caps the
                transfer.
            key: path to the Aspera private key; discovered if omitted.
        """
        if not ascp_available():
            raise BackendError(
                "ascp was not found on PATH. Install the Aspera client with: "
                "conda install -c hcc aspera-cli -- or use --backend aria2."
            )
        self.key = key or find_key()
        if not self.key:
            raise BackendError(
                "ascp is installed but its private key could not be found. "
                "Set ASPERA_SCP_KEY to the path of asperaweb_id_dsa.openssh."
            )
        self.target_rate = target_rate
        self.port = port
        self.quiet = quiet

    def build_command(self, file: RemoteFile, dest_dir: Path) -> list[str]:
        """Return the ascp argument vector for one file."""
        target = Path(file.target)
        directory = dest_dir / target.parent if str(target.parent) != "." else dest_dir
        return [
            ASCP,
            "-QT",              # fair-transfer policy, no encryption of payload
            "-l", self.target_rate,
            "-P", str(self.port),
            "-k", "1",          # resume by file attributes
            "-i", self.key,
            remote_spec(file),
            str(directory / target.name),
        ]

    def download(self, files: Sequence[RemoteFile], dest_dir: Path) -> Path:
        """Download *files* sequentially with ascp.

        Aspera saturates a link on its own, so files are transferred one at a
        time rather than in parallel.
        """
        dest_dir = Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)

        missing = [f.filename for f in files if not f.aspera]
        if missing:
            raise BackendError(
                f"{len(missing)} file(s) have no Aspera endpoint "
                f"(first: {missing[0]}). Use --backend aria2 for these."
            )

        for f in files:
            target = Path(f.target)
            (dest_dir / target.parent).mkdir(parents=True, exist_ok=True)
            proc = subprocess.run(
                self.build_command(f, dest_dir),
                capture_output=self.quiet,
                text=True,
            )
            if proc.returncode != 0:
                detail = (proc.stderr or "").strip()[:300] if self.quiet else ""
                raise BackendError(
                    f"ascp failed for {f.filename} (exit {proc.returncode}). "
                    f"If this hangs or times out, outbound UDP 33001 is "
                    f"probably blocked by a firewall. {detail}"
                )
        return dest_dir
