"""Pure-Python transfer backend.

Used when aria2c is unavailable. It downloads several files at once and
resumes partial files with a Range request, but does not split a single file
across connections -- for large files aria2c is substantially faster.
"""

from __future__ import annotations

import hashlib
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

from ..errors import BackendError, ChecksumError
from ..http import make_session
from ..models import RemoteFile

CHUNK = 1 << 20


def file_md5(path: Path) -> str:
    """Return the md5 of *path*, read in chunks so large files stay out of RAM."""
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


class RequestsBackend:
    """Threaded downloader built on :mod:`requests`."""

    name = "python"

    def __init__(
        self,
        *,
        concurrent: int = 4,
        max_tries: int = 10,
        quiet: bool = False,
        session: requests.Session | None = None,
        **_: object,
    ) -> None:
        self.concurrent = max(1, concurrent)
        # Multi-gigabyte transfers get dropped mid-stream often enough that a
        # handful of attempts is not enough; each retry resumes from the
        # partial file, so attempts are cheap and make forward progress.
        self.max_tries = max_tries
        self.quiet = quiet
        self.session = session or make_session()

    def _download_one(self, f: RemoteFile, dest_dir: Path) -> Path:
        target = dest_dir / f.target
        target.parent.mkdir(parents=True, exist_ok=True)
        part = target.with_suffix(target.suffix + ".part")

        if target.exists() and f.size and target.stat().st_size == f.size:
            return target

        offset = part.stat().st_size if part.exists() else 0
        headers = {"Range": f"bytes={offset}-"} if offset else {}

        response = self.session.get(
            f.url, headers=headers, stream=True, timeout=120
        )
        if offset and response.status_code == 200:
            # Server ignored the Range header; start over.
            offset = 0
            part.unlink(missing_ok=True)
        response.raise_for_status()

        mode = "ab" if offset else "wb"
        with open(part, mode) as handle:
            for chunk in response.iter_content(CHUNK):
                if chunk:
                    handle.write(chunk)

        if f.size and part.stat().st_size != f.size:
            raise BackendError(
                f"{f.filename}: expected {f.size} bytes, got {part.stat().st_size}."
            )
        if f.md5:
            actual = file_md5(part)
            if actual != f.md5.lower():
                part.unlink(missing_ok=True)
                raise ChecksumError(
                    f"{f.filename}: md5 mismatch (expected {f.md5}, got {actual}). "
                    f"The partial file was removed so the next run starts clean."
                )
        part.replace(target)
        return target

    def download(self, files: Sequence[RemoteFile], dest_dir: Path) -> Path:
        """Download *files* into *dest_dir*, several at a time."""
        dest_dir = Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)
        if not files:
            return dest_dir

        errors: list[str] = []
        with ThreadPoolExecutor(max_workers=self.concurrent) as pool:
            futures = {pool.submit(self._attempt, f, dest_dir): f for f in files}
            for future in as_completed(futures):
                f = futures[future]
                try:
                    future.result()
                except Exception as exc:
                    errors.append(f"{f.filename}: {exc}")
                else:
                    if not self.quiet:
                        print(f"  done  {f.target}")

        if errors:
            raise BackendError(
                f"{len(errors)} of {len(files)} downloads failed:\n  "
                + "\n  ".join(errors)
            )
        return dest_dir

    def _attempt(self, f: RemoteFile, dest_dir: Path) -> Path:
        last: Exception | None = None
        for _ in range(self.max_tries):
            try:
                return self._download_one(f, dest_dir)
            except ChecksumError:
                # A bad checksum is a corrupt transfer, not a flaky one; the
                # server will return the same bytes again.
                raise
            except Exception as exc:
                last = exc
        raise BackendError(str(last))
