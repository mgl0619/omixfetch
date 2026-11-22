"""Deciding from a BAM header whether it converts back to complete FASTQ."""

import pytest

from omixfetch.models import RemoteFile
from omixfetch.pipeline import RunPlan, StateStore, Tools, process_run
from omixfetch.probe import (
    COMPLETE,
    IMPOSSIBLE,
    MAPPED_ONLY,
    UNKNOWN,
    BamAssessment,
    assess_bam,
)

TENX = "@HD\tVN:1.4\n@CO\t10x_bam_to_fastq:R1(CR:CY,UR:UY)\n@CO\t10x_bam_to_fastq:R2(SEQ:QUAL)\n"


def pg(name, cl):
    return f"@HD\tVN:1.6\n@PG\tID:{name}\tPN:{name}\tVN:1\tCL:{cl}\n"


@pytest.mark.parametrize(
    "header,converter,complete",
    [
        # GSE184708's BAM: Cell Ranger, with explicit rebuild instructions.
        (TENX, "bamtofastq", COMPLETE),
        (pg("STAR", "STAR --genomeDir g --outSAMunmapped Within --readFilesIn a b"),
         "samtools", COMPLETE),
        (pg("STAR", "STAR --genomeDir g --readFilesIn a b"), "samtools", MAPPED_ONLY),
        (pg("hisat2", "hisat2 -x idx -1 a -2 b"), "samtools", COMPLETE),
        (pg("hisat2", "hisat2 -x idx -1 a -2 b --no-unal"), "samtools", MAPPED_ONLY),
        (pg("bowtie2", "bowtie2 -x idx -U a --no-unal"), "samtools", MAPPED_ONLY),
        (pg("bwa", "bwa mem ref a b"), "samtools", COMPLETE),
        (pg("minimap2", "minimap2 -a ref a"), "samtools", COMPLETE),
        (pg("minimap2", "minimap2 -a --sam-hit-only ref a"), "samtools", MAPPED_ONLY),
        (pg("mystery", "mystery --in a"), "samtools", UNKNOWN),
        # A single-cell BAM without 10x's lines: barcodes cannot become reads.
        (pg("STAR", "STAR --soloType CB_UMI_Simple --outSAMunmapped Within"), None, IMPOSSIBLE),
    ],
)
def test_assess_bam(header, converter, complete):
    verdict = assess_bam(header)
    assert verdict.converter == converter
    assert verdict.complete == complete
    assert verdict.reason


def test_only_complete_conversions_are_acceptable():
    assert BamAssessment("samtools", COMPLETE, "bwa", "").acceptable
    assert not BamAssessment("samtools", MAPPED_ONLY, "star", "").acceptable
    assert not BamAssessment("samtools", UNKNOWN, None, "").acceptable
    assert not BamAssessment(None, IMPOSSIBLE, None, "").acceptable


# ---- pipeline behaviour ---------------------------------------------------

def _plan():
    bam = RemoteFile(url="https://sra-pub-src-2.s3.amazonaws.com/SRR1/a.bam",
                     filename="a.bam", size=10)
    return RunPlan("SRR1", "GSM1", None, {"original": [bam]})


def _tools(verdict, fetched, converters=frozenset({"bamtofastq", "samtools"})):
    def download(files, dest):
        fetched.append(files[0].filename)
        dest.mkdir(parents=True, exist_ok=True)
        for f in files:
            (dest / f.filename).write_bytes(b"x")

    used = {}

    def convert(bam, out, converter):
        used["converter"] = converter
        fq = out / "a_1.fastq.gz"
        fq.write_bytes(b"x")
        return [fq]

    tools = Tools(download=download, convert_bam=convert, check_bam=lambda url: verdict,
                  converters=converters, log=lambda m: None)
    return tools, used


def _run(tmp_path, tools, **kw):
    return process_run(_plan(), tmp_path, StateStore(tmp_path / "s.json"), tools,
                       check_space=False, **kw)


def test_complete_bam_uses_the_assessed_converter(tmp_path):
    fetched = []
    tools, used = _tools(BamAssessment("samtools", COMPLETE, "bwa", "ok"), fetched)
    assert _run(tmp_path, tools).status == "done"
    assert used["converter"] == "samtools"


def test_mapped_only_bam_is_refused_before_download(tmp_path):
    fetched = []
    tools, _ = _tools(BamAssessment("samtools", MAPPED_ONLY, "star", "drops unmapped"), fetched)
    state = _run(tmp_path, tools)
    assert state.status == "failed" and "mapped-only" in state.error
    assert fetched == []


def test_mapped_only_bam_accepted_when_allowed(tmp_path):
    fetched = []
    tools, _ = _tools(BamAssessment("samtools", MAPPED_ONLY, "star", "drops unmapped"), fetched)
    assert _run(tmp_path, tools, allow_incomplete=True).status == "done"


def test_impossible_bam_is_refused_even_when_allowed(tmp_path):
    fetched = []
    tools, _ = _tools(BamAssessment(None, IMPOSSIBLE, "star", "no barcodes"), fetched)
    state = _run(tmp_path, tools, allow_incomplete=True)
    assert state.status == "failed" and fetched == []


def test_missing_converter_is_refused_before_download(tmp_path):
    fetched = []
    tools, _ = _tools(BamAssessment("samtools", COMPLETE, "bwa", "ok"), fetched,
                      converters=frozenset({"bamtofastq"}))
    state = _run(tmp_path, tools)
    assert state.status == "failed" and "samtools" in state.error
    assert fetched == []
