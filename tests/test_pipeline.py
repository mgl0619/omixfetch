"""Autonomous GEO -> FASTQ pipeline, exercised with fake tools offline."""

from pathlib import Path
from types import SimpleNamespace

import pytest

from omixfetch.errors import BackendError
from omixfetch.models import RemoteFile
from omixfetch.pipeline import (
    DONE,
    FAILED,
    SKIPPED,
    RunPlan,
    StateStore,
    Tools,
    check_disk,
    disk_needed,
    is_complete,
    process_run,
    verify_spots,
    write_manifest,
)


def _file(name, size=10):
    return RemoteFile(url=f"https://x/{name}", filename=name, size=size)


def _writer(calls):
    """A download that creates each file and records what was fetched."""
    def download(files, dest):
        calls.append([f.filename for f in files])
        dest.mkdir(parents=True, exist_ok=True)
        for f in files:
            (dest / f.filename).write_bytes(b"x")
    return download


def _sra_converter(spots):
    def convert(sra, out):
        fq = out / f"{sra.stem}_1.fastq.gz"
        fq.write_bytes(b"@r\nA\n+\nI\n")
        return SimpleNamespace(files=[fq], spots_read=spots)
    return convert


@pytest.fixture
def store(tmp_path):
    return StateStore(tmp_path / "state.json")


def test_ena_fastq_route_needs_no_conversion(tmp_path, store):
    calls = []
    plan = RunPlan("SRR1", "GSM1", 100, {"ena-fastq": [_file("SRR1_1.fastq.gz")]})
    state = process_run(plan, tmp_path, store, Tools(download=_writer(calls), log=lambda m: None),
                        check_space=False)
    assert state.status == DONE and state.route == "ena-fastq"
    assert Path(state.fastq[0]) == tmp_path / "fastq" / "GSM1" / "SRR1_1.fastq.gz"
    assert Path(state.fastq[0]).exists()


def test_sra_route_converts_verifies_and_cleans_up(tmp_path, store):
    plan = RunPlan("SRR2", "GSM2", 500, {"sra": [_file("SRR2.sra")]})
    tools = Tools(download=_writer([]), convert_sra=_sra_converter(500), log=lambda m: None)
    state = process_run(plan, tmp_path, store, tools, check_space=False)
    assert state.status == DONE and state.spots_written == 500
    assert (tmp_path / "fastq" / "GSM2" / "SRR2_1.fastq.gz").exists()
    # The .sra and its work directory are removed after a verified conversion.
    assert not (tmp_path / "work" / "SRR2").exists()


def test_falls_through_to_next_route_on_failure(tmp_path, store):
    calls = []

    def flaky(files, dest):
        if files[0].filename.endswith(".sra"):
            raise BackendError("reference fetch failed")
        _writer(calls)(files, dest)

    plan = RunPlan("SRR3", "GSM3", None, {
        "sra": [_file("SRR3.sra")],
        "original": [_file("SRR3.bam")],
    })
    def bam(bam_path, out):
        fq = out / "SRR3_R1.fastq.gz"
        fq.write_bytes(b"x")
        return [fq]
    tools = Tools(download=flaky, convert_sra=_sra_converter(1), convert_bam=bam,
                  log=lambda m: None)
    state = process_run(plan, tmp_path, store, tools, check_space=False)
    assert state.status == DONE and state.route == "original"
    assert state.attempts == ["sra", "original"]


def test_spot_mismatch_rejects_the_route(tmp_path, store):
    """A truncated conversion must not be accepted as finished."""
    plan = RunPlan("SRR4", "GSM4", 1000, {"sra": [_file("SRR4.sra")]})
    tools = Tools(download=_writer([]), convert_sra=_sra_converter(999), log=lambda m: None)
    state = process_run(plan, tmp_path, store, tools, check_space=False)
    assert state.status == FAILED
    assert "spot count mismatch" in state.error


def test_missing_converter_is_refused_before_downloading(tmp_path, store):
    """Never spend a large transfer on a route that cannot finish."""
    calls = []
    plan = RunPlan("SRR5", "GSM5", None, {"sra": [_file("SRR5.sra", 10**11)]})
    state = process_run(plan, tmp_path, store,
                        Tools(download=_writer(calls), convert_sra=None, log=lambda m: None),
                        check_space=False)
    assert state.status == FAILED
    assert calls == []


def test_completed_runs_are_skipped_on_rerun(tmp_path, store):
    calls = []
    plan = RunPlan("SRR6", "GSM6", None, {"ena-fastq": [_file("SRR6.fastq.gz")]})
    tools = Tools(download=_writer(calls), log=lambda m: None)
    process_run(plan, tmp_path, store, tools, check_space=False)
    process_run(plan, tmp_path, StateStore(store.path), tools, check_space=False)
    assert len(calls) == 1


def test_done_run_is_refetched_if_its_files_vanished(tmp_path, store):
    calls = []
    plan = RunPlan("SRR7", "GSM7", None, {"ena-fastq": [_file("SRR7.fastq.gz")]})
    tools = Tools(download=_writer(calls), log=lambda m: None)
    state = process_run(plan, tmp_path, store, tools, check_space=False)
    Path(state.fastq[0]).unlink()
    assert not is_complete(store.get("SRR7"))
    process_run(plan, tmp_path, store, tools, check_space=False)
    assert len(calls) == 2


def test_run_with_no_public_source_is_skipped(tmp_path, store):
    state = process_run(RunPlan("SRR8", "GSM8", None, {}), tmp_path, store,
                        Tools(download=_writer([]), log=lambda m: None), check_space=False)
    assert state.status == SKIPPED
    assert "controlled access" in state.error


def test_state_survives_a_restart(tmp_path):
    first = StateStore(tmp_path / "s.json")
    first.update("SRR9", status=DONE, fastq=["/a"])
    assert StateStore(tmp_path / "s.json").get("SRR9").status == DONE


def test_verify_spots():
    assert verify_spots(10, 10) is None
    assert verify_spots(None, 10) is None
    assert verify_spots(10, None) is None
    assert "mismatch" in verify_spots(10, 9)


def test_disk_estimate_is_largest_for_sra():
    plan = RunPlan("R", None, None, {"sra": [_file("a", 100)], "ena-fastq": [_file("b", 100)]})
    assert disk_needed(plan, "sra") > disk_needed(plan, "ena-fastq")


def test_check_disk_refuses_impossible_requests(tmp_path):
    with pytest.raises(BackendError, match="not enough disk"):
        check_disk(tmp_path, 10**18)


def test_manifest_lists_every_run(tmp_path, store):
    store.update("SRR1", sample="GSM1", status=DONE, route="sra",
                 fastq=["/f/a_1.fastq.gz", "/f/a_2.fastq.gz"])
    store.update("SRR2", sample="GSM2", status=FAILED, error="boom")
    lines = write_manifest(store, tmp_path / "m.tsv").read_text().splitlines()
    assert lines[0].startswith("sample\trun\tstatus")
    assert len(lines) == 3
    assert "a_1.fastq.gz;/f/a_2.fastq.gz" in lines[1]
