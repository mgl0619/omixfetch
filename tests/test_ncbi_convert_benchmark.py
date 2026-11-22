"""NCBI SDL resolution, conversion commands and throughput sampling."""

import pytest

from omixfetch.benchmark import SampleResult, benchmark_source, format_table
from omixfetch.convert import PATH_ENV_VAR, build_command, resolve_fasterq_dump
from omixfetch.models import RemoteFile
from omixfetch.sources.ncbi import objects_to_files, parse_sdl

SDL = {
    "version": "2",
    "result": [
        {
            "bundle": "SRR16036138",
            "status": 200,
            "msg": "ok",
            "files": [
                {
                    "object": "remote|127955708",
                    "accession": "SRR16036138",
                    "type": "TenX",
                    "name": "170830_AGP_XX_E10.bam",
                    "size": 119194087864,
                    "md5": "00dbf15405d748340fbbf71240bf6771",
                    "locations": [
                        {"service": "s3", "region": "us-east-1",
                         "link": "https://sra-pub-src-2.s3.amazonaws.com/SRR16036138/x.bam.1"}
                    ],
                },
                {
                    "object": "srapub|SRR16036138",
                    "accession": "SRR16036138",
                    "type": "sra",
                    "name": "SRR16036138",
                    "size": 48085164879,
                    "md5": "b5a9d3e109005cd980d5a1d4e32abbf2",
                    "locations": [
                        {"service": "s3", "region": "us-east-1",
                         "link": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR16036138/SRR16036138"}
                    ],
                },
            ],
        },
        # A restricted run must not abort the whole lookup.
        {"bundle": "SRR000000", "status": 403, "msg": "forbidden", "files": []},
    ],
}


def test_parse_sdl_classifies_objects_and_skips_failures():
    objects = parse_sdl(SDL)
    assert len(objects) == 2
    assert {o.flavour for o in objects} == {"sra", "original"}

    sra = next(o for o in objects if o.flavour == "sra")
    assert sra.size == 48085164879
    assert sra.md5 == "b5a9d3e109005cd980d5a1d4e32abbf2"

    original = next(o for o in objects if o.flavour == "original")
    # The .sra is substantially smaller than the submitted BAM.
    assert original.size > sra.size * 2


def test_parse_sdl_ignores_records_without_a_location():
    payload = {"result": [{"bundle": "X", "status": 200,
                           "files": [{"object": "srapub|X", "name": "X", "locations": []}]}]}
    assert parse_sdl(payload) == []


def test_objects_to_files_selects_flavour_and_names_sra():
    files = objects_to_files(parse_sdl(SDL), flavour="sra")
    assert len(files) == 1
    # srapub objects are named for the run; an extension is added so
    # downstream tools recognise the file.
    assert files[0].filename == "SRR16036138.sra"
    assert files[0].kind == "raw"
    assert files[0].md5


def test_objects_to_files_maps_runs_to_samples():
    files = objects_to_files(parse_sdl(SDL), flavour="original",
                             sample_by_run={"SRR16036138": "GSM1"})
    assert files[0].accession == "GSM1"


def test_objects_to_files_rejects_unknown_flavour():
    with pytest.raises(ValueError, match="Unknown flavour"):
        objects_to_files(parse_sdl(SDL), flavour="bam")


def test_fasterq_dump_command(tmp_path):
    cmd = build_command(tmp_path / "SRR1.sra", tmp_path / "out", threads=8)
    assert cmd[0] == "fasterq-dump"
    assert "--split-3" in cmd
    assert "--threads" in cmd and "8" in cmd
    # Scratch goes beside the output, not /tmp, which is small on many nodes.
    assert "--temp" in cmd


def test_fasterq_dump_rejects_unknown_split(tmp_path):
    with pytest.raises(ValueError, match="Unknown split mode"):
        build_command(tmp_path / "a.sra", tmp_path, split="pairs")


def test_benchmark_reports_rate_and_projection():
    f = RemoteFile(url="https://x/a", filename="a")
    result = benchmark_source("src", f, session=None,
                              sampler=lambda s, u, n: (100 * 1024**2, 2.0))
    assert result.mib_per_s == pytest.approx(50.0)
    # One hour at 50 MiB/s is exactly 50 * 3600 MiB.
    assert result.projected_hours(50 * 3600 * 1024**2) == pytest.approx(1.0)


def test_benchmark_turns_failures_into_results():
    f = RemoteFile(url="https://x/a", filename="a")

    def boom(*a, **k):
        raise ConnectionError("blocked")

    result = benchmark_source("src", f, session=None, sampler=boom)
    assert result.error and result.mib_per_s is None


def test_format_table_ranks_fastest_first_and_lists_failures_last():
    rows = format_table(
        [
            SampleResult("slow", "a", 1024**2, 1.0),
            SampleResult("broken", "b", 0, 0.0, error="nope"),
            SampleResult("fast", "c", 100 * 1024**2, 1.0),
        ],
        None,
    )
    assert rows[0].startswith("fast")
    assert rows[1].startswith("slow")
    assert rows[2].startswith("broken")


def _fake_binary(tmp_path, name="fasterq-dump"):
    path = tmp_path / name
    path.write_text("#!/bin/sh\n")
    path.chmod(0o755)
    return path


def test_resolve_prefers_explicit_path(tmp_path):
    binary = _fake_binary(tmp_path)
    assert resolve_fasterq_dump(str(binary)) == str(binary)


def test_resolve_accepts_a_directory(tmp_path):
    binary = _fake_binary(tmp_path)
    assert resolve_fasterq_dump(str(tmp_path)) == str(binary)


def test_resolve_falls_back_to_environment_variable(tmp_path, monkeypatch):
    """The toolkit need not live in omixfetch's own environment."""
    binary = _fake_binary(tmp_path)
    monkeypatch.setenv(PATH_ENV_VAR, str(binary))
    assert resolve_fasterq_dump(None) == str(binary)


def test_resolve_ignores_a_nonexistent_path(tmp_path, monkeypatch):
    monkeypatch.delenv(PATH_ENV_VAR, raising=False)
    monkeypatch.setattr("shutil.which", lambda _n: None)
    assert resolve_fasterq_dump(str(tmp_path / "absent")) is None


def test_build_command_uses_the_resolved_binary(tmp_path):
    cmd = build_command(tmp_path / "a.sra", tmp_path, binary="/opt/sra/fasterq-dump")
    assert cmd[0] == "/opt/sra/fasterq-dump"
