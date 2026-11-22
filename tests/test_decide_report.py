"""The route decision tree and the HTML report."""

from omixfetch.decide import BAM, FASTQ, decide, submission_format
from omixfetch.models import RemoteFile
from omixfetch.pipeline import DONE, FAILED, RunPlan, StateStore
from omixfetch.probe import COMPLETE, IMPOSSIBLE, MAPPED_ONLY, BamAssessment
from omixfetch.report import write_report

GIB = 1024**3


def _f(name, size):
    return RemoteFile(url=f"https://x/{name}", filename=name, size=size)


# ---- what was submitted -----------------------------------------------------

def test_submission_format_from_ncbi_original_files():
    assert submission_format(["170830_AGP_XX_E10.bam"], "BAM;BAI") == BAM
    assert submission_format(["a.cram"], None) == "cram"
    # FASTQ submissions have no original file listed at NCBI.
    assert submission_format([], "") == FASTQ
    assert submission_format([], "FASTQ") == FASTQ
    # ENA as a fallback when NCBI lists nothing.
    assert submission_format([], "BAM") == BAM


# ---- BAM submissions --------------------------------------------------------

def _bam_plan(verdict, mates_missing=False):
    return RunPlan(
        "SRR16036138", "GSM5595982", 1_091_833_030,
        {"sra": [_f("SRR16036138.sra", 45 * GIB)],
         "original": [_f("170830_AGP_XX_E10.bam", 111 * GIB)]},
        mates_missing=mates_missing, submitted=BAM, bam_verdict=verdict,
    )


def test_tenx_bam_uses_bamtofastq_and_excludes_derived_routes():
    """GSE184708: the smaller .sra is excluded because it lost R1."""
    d = decide(_bam_plan(BamAssessment("bamtofastq", COMPLETE, "cellranger", "10x"),
                         mates_missing=True))
    assert d.order == ["original"]
    assert "barcode" in d.excluded["sra"]
    assert "10x" in d.branch


def test_complete_generic_bam_prefers_the_bam():
    d = decide(_bam_plan(BamAssessment("samtools", COMPLETE, "bwa", "keeps unmapped")))
    assert d.order[0] == "original"


def test_lossy_bam_has_no_faithful_route():
    d = decide(_bam_plan(BamAssessment("samtools", MAPPED_ONLY, "star", "drops unmapped")))
    assert d.order == []
    assert "no faithful route" in d.branch
    assert set(d.excluded) >= {"original", "sra"}


def test_lossy_routes_only_when_explicitly_allowed():
    d = decide(_bam_plan(BamAssessment(None, IMPOSSIBLE, "star", "no barcodes")),
               allow_incomplete=True)
    assert d.lossy and d.order  # something to try, flagged as lossy


# ---- FASTQ submissions ------------------------------------------------------

def _fastq_plan():
    return RunPlan(
        "SRR19429286", "GSM6203352", 23959,
        {"ena-fastq": [_f("a_1.fastq.gz", 800_000), _f("a_2.fastq.gz", 800_000)],
         "sra": [_f("SRR19429286.sra", 1_600_000)]},
        submitted=FASTQ,
    )


def test_fastq_submission_picks_the_faster_source():
    # Your case: ENA slow, NCBI fast -> the .sra wins despite conversion.
    d = decide(_fastq_plan(), rates={"ena": 8.0, "ncbi": 55.0})
    assert d.order[0] == "sra"
    assert d.estimates["sra"] < d.estimates["ena-fastq"]


def test_fastq_submission_prefers_no_conversion_when_rates_are_equal():
    d = decide(_fastq_plan(), rates={"ena": 50.0, "ncbi": 50.0})
    assert d.order[0] == "ena-fastq"


def test_user_order_overrides_ranking_but_not_exclusions():
    d = decide(_bam_plan(BamAssessment("bamtofastq", COMPLETE, "cellranger", "10x")),
               order_override=("sra", "original"))
    assert d.order == ["original"]  # sra stays excluded: it lost R1


def test_no_public_files():
    d = decide(RunPlan("SRR0", "GSM0", None, {}))
    assert d.order == [] and "controlled access" in d.branch


# ---- report ----------------------------------------------------------------

def test_report_is_self_contained_and_lists_every_run(tmp_path):
    plans = [_bam_plan(BamAssessment("bamtofastq", COMPLETE, "cellranger", "10x")),
             _fastq_plan()]
    for plan in plans:
        plan.decision = decide(plan, rates={"ena": 8.0, "ncbi": 55.0})
    store = StateStore(tmp_path / "s.json")
    store.update("SRR16036138", status=DONE, route="original", bytes=111 * GIB,
                 started=0.0, finished=3600.0, fastq=["/o/a_R1_001.fastq.gz"])
    store.update("SRR19429286", status=FAILED, attempts=["sra"], error="boom <b>")

    html = write_report(tmp_path / "r.html", "GSE184708", plans, store,
                        title="Gonad development", rates={"ena": 8.0, "ncbi": 55.0},
                        tools={"bamtofastq": "/t/bamtofastq"}).read_text()
    assert "SRR16036138" in html and "SRR19429286" in html
    assert "lost the barcode" in html          # the exclusion reason is shown
    assert "boom &lt;b&gt;" in html            # errors are escaped
    assert "<script" not in html and "http" not in html.split("<body>")[1].split("Tools")[0]
    assert "55.0 MiB/s" in html
