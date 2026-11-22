"""Live tests against the real archives, in cost tiers.

Skipped unless OMIXFETCH_LIVE names a tier; each tier includes the cheaper ones.

    OMIXFETCH_LIVE=quick   python -m pytest tests/live    # ~1 min, ~10 MB
    OMIXFETCH_LIVE=medium  python -m pytest tests/live    # ~30 min, ~25 GB disk
    OMIXFETCH_LIVE=final   python -m pytest tests/live    # many hours, 1.7 TiB

GSE184708 belongs to the final tier only: never iterate against it.
"""

import gzip
import os
import subprocess
import sys
from pathlib import Path

import pytest

TIERS = ("quick", "medium", "final")
LEVEL = os.environ.get("OMIXFETCH_LIVE", "")


def tier(name):
    enabled = LEVEL in TIERS and TIERS.index(LEVEL) >= TIERS.index(name)
    return pytest.mark.skipif(not enabled, reason=f"set OMIXFETCH_LIVE={name} to run")


def omixfetch(*args, cwd=None):
    return subprocess.run(
        [sys.executable, "-m", "omixfetch", *args],
        capture_output=True, text=True, cwd=cwd,
    )


def manifest(path):
    lines = Path(path).read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=False)) for line in lines[1:]]


def reads(path):
    with gzip.open(path, "rt") as handle:
        return sum(1 for _ in handle) // 4


# Two runs of a few MB each from GSE204978; ENA serves FASTQ for both and
# their .sra archives are unaligned, so conversion needs no reference.
SMALL = ("GSE204978", ["SRR19429286", "SRR19429291"], {"SRR19429286": 23959,
                                                        "SRR19429291": 26554})


def _include(runs):
    return [arg for run in runs for arg in ("--include", run)]


@tier("quick")
@pytest.mark.parametrize("route", ["ena-fastq", "sra"])
def test_small_series_each_route(tmp_path, route):
    gse, runs, spots = SMALL
    result = omixfetch("fastq", gse, "-o", str(tmp_path), "--routes", route, *_include(runs))
    assert result.returncode == 0, result.stdout + result.stderr

    rows = {r["run"]: r for r in manifest(tmp_path / gse / "fastq_manifest.tsv")}
    for run in runs:
        assert rows[run]["status"] == "done" and rows[run]["route"] == route
        files = rows[run]["fastq"].split(";")
        assert len(files) == 2  # paired
        for f in files:
            assert reads(f) == spots[run]
    assert not (tmp_path / gse / "work" / runs[0]).exists()


@tier("quick")
def test_rerun_skips_finished_runs(tmp_path):
    gse, runs, _ = SMALL
    omixfetch("fastq", gse, "-o", str(tmp_path), *_include(runs))
    again = omixfetch("fastq", gse, "-o", str(tmp_path), *_include(runs))
    assert again.returncode == 0
    assert (again.stdout + again.stderr).count("already done, skipping") == len(runs)


@tier("quick")
def test_controlled_access_series_is_reported(tmp_path):
    result = omixfetch("fastq", "GSE164378", "-o", str(tmp_path), "--dry-run")
    assert result.returncode == 1
    assert "no public runs" in result.stdout + result.stderr


@tier("medium")
def test_tenx_bam_converts_with_barcodes(tmp_path):
    """Smallest 10x BAM found: GSE220167 / SRR22543607, 8.3 GB."""
    gse, run = "GSE220167", "SRR22543607"
    result = omixfetch("fastq", gse, "-o", str(tmp_path), "--include", run, "-t", "8")
    assert result.returncode == 0, result.stdout + result.stderr
    row = next(r for r in manifest(tmp_path / gse / "fastq_manifest.tsv") if r["run"] == run)
    assert row["route"] == "original"  # the .sra would have lost R1
    names = " ".join(Path(f).name for f in row["fastq"].split(";"))
    for read in ("_R1_", "_R2_"):
        assert read in names, names
    r1 = next(f for f in row["fastq"].split(";") if "_R1_" in Path(f).name)
    with gzip.open(r1, "rt") as handle:
        handle.readline()
        barcode_umi = handle.readline().strip()
    assert len(barcode_umi) >= 26  # 16 bp barcode + 10-12 bp UMI


@tier("final")
def test_gse184708_full(tmp_path):
    out = Path(os.environ.get("OMIXFETCH_LIVE_OUT", tmp_path))
    result = omixfetch("fastq", "GSE184708", "-o", str(out), "-t", "16")
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]
