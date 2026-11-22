"""Measure achievable throughput from each available source.

Which mirror is fastest depends entirely on the network path of the machine
doing the downloading -- a cluster node, a laptop and a cloud VM can differ
by an order of magnitude, and published benchmarks are therefore useless.
This module samples a fixed number of bytes from one representative file per
source and reports the rate, so the decision can be made with local numbers.

Sampling uses HTTP range requests, so it costs the sample size rather than
the size of the file.
"""

from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import dataclass

import requests

from .models import RemoteFile

#: Bytes sampled per source unless overridden.
DEFAULT_SAMPLE = 100 * 1024 * 1024

CHUNK = 1 << 20


@dataclass
class SampleResult:
    """Throughput measured for one source."""

    source: str
    filename: str
    bytes_read: int
    seconds: float
    error: str | None = None

    @property
    def mib_per_s(self) -> float | None:
        if self.error or self.seconds <= 0:
            return None
        return self.bytes_read / self.seconds / 1024**2

    def projected_hours(self, total_bytes: int) -> float | None:
        """Hours to move *total_bytes* at the measured rate."""
        rate = self.mib_per_s
        if not rate:
            return None
        return total_bytes / (rate * 1024**2) / 3600


def sample_http(
    session: requests.Session,
    url: str,
    sample_bytes: int = DEFAULT_SAMPLE,
    timeout: int = 180,
) -> tuple[int, float]:
    """Read up to *sample_bytes* from *url* and return ``(bytes, seconds)``.

    The clock starts after the response headers arrive, so DNS and TLS setup
    are excluded and the figure reflects sustained transfer rate.
    """
    headers = {"Range": f"bytes=0-{sample_bytes - 1}"}
    response = session.get(url, headers=headers, stream=True, timeout=timeout)
    response.raise_for_status()

    read = 0
    start = time.perf_counter()
    for chunk in response.iter_content(CHUNK):
        read += len(chunk)
        if read >= sample_bytes:
            break
    elapsed = time.perf_counter() - start
    response.close()
    return read, elapsed


def benchmark_source(
    name: str,
    file: RemoteFile,
    session: requests.Session,
    sample_bytes: int = DEFAULT_SAMPLE,
    sampler: Callable[..., tuple[int, float]] | None = None,
) -> SampleResult:
    """Measure one source, converting failures into a reported error."""
    sampler = sampler or sample_http
    try:
        read, elapsed = sampler(session, file.url, sample_bytes)
    except Exception as exc:  # network failures are a result, not a crash
        return SampleResult(name, file.filename, 0, 0.0, error=str(exc)[:160])
    return SampleResult(name, file.filename, read, elapsed)


def format_table(results: list[SampleResult], total_bytes: int | None) -> list[str]:
    """Render results as aligned text rows, fastest first."""
    ranked = sorted(
        results,
        key=lambda r: (r.mib_per_s is None, -(r.mib_per_s or 0)),
    )
    rows = []
    for r in ranked:
        if r.error:
            rows.append(f"{r.source:<22} failed: {r.error[:60]}")
            continue
        line = (
            f"{r.source:<22} {r.mib_per_s:>7.1f} MiB/s  "
            f"({r.bytes_read / 1e6:.0f} MB in {r.seconds:.1f}s)"
        )
        if total_bytes:
            hours = r.projected_hours(total_bytes)
            line += f"  -> {hours:>6.1f} h for the full set"
        rows.append(line)
    return rows
