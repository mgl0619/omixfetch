"""Bundled tool installation, BAM header probing and RunInfo parsing."""

import gzip
import io
import tarfile

import pytest

from omixfetch import bundled
from omixfetch.errors import BackendError
from omixfetch.models import RemoteFile
from omixfetch.pipeline import RunPlan, StateStore, Tools, process_run, verify_spots
from omixfetch.probe import bam_header, is_tenx_bam
from omixfetch.sources.runinfo import parse_runinfo

# ---- platform selection -------------------------------------------------

@pytest.mark.parametrize(
    "system,machine,build",
    [("Darwin", "arm64", "mac-arm64"), ("Darwin", "x86_64", "mac64")],
)
def test_sra_toolkit_build_for_macs(system, machine, build):
    assert bundled.sra_toolkit_build(system, machine) == build


def test_sra_toolkit_build_for_linux(monkeypatch):
    monkeypatch.setattr(bundled, "_linux_flavour", lambda: "alma_linux64")
    assert bundled.sra_toolkit_build("Linux", "x86_64") == "alma_linux64"


def test_sra_toolkit_build_rejects_unsupported():
    with pytest.raises(BackendError, match="no prebuilt"):
        bundled.sra_toolkit_build("Linux", "aarch64")


def test_bamtofastq_asset():
    assert bundled.bamtofastq_asset("Darwin") == "bamtofastq_macos"
    assert bundled.bamtofastq_asset("Linux") == "bamtofastq_linux"


def test_home_is_overridable(tmp_path, monkeypatch):
    monkeypatch.setenv(bundled.HOME_VAR, str(tmp_path))
    assert bundled.bin_dir() == tmp_path / "tools" / "bin"


def test_bundled_finds_an_installed_tool(tmp_path, monkeypatch):
    monkeypatch.setenv(bundled.HOME_VAR, str(tmp_path))
    assert bundled.bundled("fasterq-dump") is None
    real = tmp_path / "real-fasterq-dump"
    real.write_text("#!/bin/sh\n")
    real.chmod(0o755)
    bundled._link(real, "fasterq-dump")
    assert bundled.bundled("fasterq-dump") == str(bundled.bin_dir() / "fasterq-dump")


def test_safe_extract_rejects_path_traversal(tmp_path):
    archive = tmp_path / "evil.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        data = b"x"
        info = tarfile.TarInfo("../escape.txt")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    with pytest.raises(BackendError, match="unsafe path"):
        bundled._safe_extract(archive, tmp_path / "dest")


# ---- BAM header probe ---------------------------------------------------

def _bam_bytes(header_text: str) -> bytes:
    """A minimal BGZF-compressed BAM containing only a header."""
    text = header_text.encode()
    payload = b"BAM\x01" + len(text).to_bytes(4, "little") + text + (0).to_bytes(4, "little")
    return gzip.compress(payload)


class _Session:
    def __init__(self, content):
        self.content = content

    def get(self, url, headers=None, timeout=None):
        return type("R", (), {"content": self.content})()


TENX = (
    "@HD\tVN:1.4\n"
    "@CO\t10x_bam_to_fastq:R1(CR:CY,UR:UY)\n"
    "@CO\t10x_bam_to_fastq:R2(SEQ:QUAL)\n"
)


def test_bam_header_reads_text_from_bgzf():
    header = bam_header(_Session(_bam_bytes(TENX)), "https://x/a.bam")
    assert "10x_bam_to_fastq:R1(CR:CY,UR:UY)" in header


def test_is_tenx_bam():
    assert is_tenx_bam(TENX)
    assert not is_tenx_bam("@HD\tVN:1.4\n@CO\tuser command line: STAR\n")


def test_bam_header_rejects_non_bam():
    with pytest.raises(ValueError, match="not a BAM"):
        bam_header(_Session(gzip.compress(b"hello world")), "https://x/a")


# ---- RunInfo ------------------------------------------------------------

RUNINFO = (
    "Run,spots,bases,avgLength,spots_with_mates,LibraryLayout\n"
    "SRR16036138,1091833030,106999636940,98,0,PAIRED\n"
    "SRR1039514,40021727,4042594427,101,40021727,PAIRED\n"
)


def test_runinfo_detects_a_lost_mate():
    info = parse_runinfo(RUNINFO)
    assert info["SRR16036138"].mates_missing  # the GSE184708 pattern
    assert not info["SRR1039514"].mates_missing
    assert info["SRR16036138"].spots == 1091833030


# ---- pipeline safeguards ------------------------------------------------

def test_zero_expected_spots_is_treated_as_unknown():
    """ENA reports 0 for BAM submissions it never processed."""
    assert verify_spots(0, 1_091_833_030) is None


def _f(name):
    return RemoteFile(url=f"https://x/{name}", filename=name, size=10)


def test_tenx_run_skips_the_lossy_sra_route(tmp_path):
    fetched = []

    def download(files, dest):
        fetched.append(files[0].filename)
        dest.mkdir(parents=True, exist_ok=True)
        for f in files:
            (dest / f.filename).write_bytes(b"x")

    def convert_bam(bam, out, converter=None):
        fq = out / "R1.fastq.gz"
        fq.write_bytes(b"x")
        return [fq]

    from omixfetch.probe import COMPLETE, BamAssessment

    plan = RunPlan("SRR1", "GSM1", None,
                   {"sra": [_f("SRR1.sra")], "original": [_f("SRR1.bam")]}, tenx=True,
                   submitted="bam",
                   bam_verdict=BamAssessment("bamtofastq", COMPLETE, "cellranger", "10x"))
    tools = Tools(download=download, convert_sra=lambda s, o: None,
                  convert_bam=convert_bam, log=lambda m: None)
    state = process_run(plan, tmp_path, StateStore(tmp_path / "s.json"), tools,
                        check_space=False)
    assert state.route == "original"
    assert fetched == ["SRR1.bam"]  # the .sra was never downloaded


def test_non_tenx_bam_is_rejected_before_download(tmp_path):
    fetched = []

    def check_bam(url):
        raise BackendError("BAM has no 10x_bam_to_fastq header lines")

    plan = RunPlan("SRR2", "GSM2", None, {"original": [_f("SRR2.bam")]})
    tools = Tools(download=lambda f, d: fetched.append(f), check_bam=check_bam,
                  convert_bam=lambda *a: [], log=lambda m: None)
    state = process_run(plan, tmp_path, StateStore(tmp_path / "s.json"), tools,
                        check_space=False)
    assert state.status == "failed"
    assert fetched == []
