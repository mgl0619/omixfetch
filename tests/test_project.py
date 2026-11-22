"""Read-archive projects resolved directly (PRJNA/SRP/SRR/SAMN ...)."""

from omixfetch.layout import apply_layout
from omixfetch.manifest import write_runs_table, write_samples_table
from omixfetch.resolve import choose_resolver, resolve
from omixfetch.sources import project
from omixfetch.sources.ena import parse_filereport

HEADER = (
    "run_accession\texperiment_accession\tsample_accession\tsecondary_sample_accession\t"
    "study_accession\tsecondary_study_accession\tsample_alias\tsample_title\t"
    "scientific_name\ttax_id\tlibrary_layout\tlibrary_strategy\tlibrary_source\t"
    "library_selection\tinstrument_model\tread_count\tbase_count\tfastq_ftp\t"
    "fastq_bytes\tfastq_md5\n"
)


def _row(run, sample, alias, md5s=("a" * 32, "b" * 32)):
    base = f"ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/0/{run}/{run}"
    return (
        f"{run}\tSRX1\t{sample}\tSRS1\tPRJNA672658\tSRP289408\t{alias}\t{alias} title\t"
        f"Homo sapiens\t9606\tPAIRED\tRNA-Seq\tTRANSCRIPTOMIC SINGLE CELL\tPCR\t"
        f"HiSeq X Ten\t100\t30000\t{base}_1.fastq.gz;{base}_2.fastq.gz\t10;20\t"
        f"{md5s[0]};{md5s[1]}\n"
    )


# Two runs from one BioSample and one from another.
FILEREPORT = (
    HEADER
    + _row("SRR12921178", "SAMN16565585", "RSA_5")
    + _row("SRR12921179", "SAMN16565584", "RSA_4")
    + _row("SRR12921180", "SAMN16565585", "RSA_5")
)

STUDY = {
    "study_accession": "PRJNA672658",
    "secondary_study_accession": "SRP289408",
    "study_title": "scRNA of normal and RSA sample",
    "study_description": "Decidua.",
    "center_name": "Xiamen University",
    "first_public": "2022-06-09",
}


def _patch(monkeypatch, runs_text=FILEREPORT, study=STUDY):
    seen = {}

    def fake_runs(session, acc, prefer="fastq"):
        seen["runs"] = acc
        return parse_filereport(runs_text, prefer=prefer)

    def fake_study(session, acc):
        seen["study"] = acc
        return dict(study)

    monkeypatch.setattr(project, "fetch_runs", fake_runs)
    monkeypatch.setattr(project, "fetch_study", fake_study)
    monkeypatch.setattr(
        project, "fetch_characteristics", lambda session, acc: {"tissue": "decidua", "acc": acc}
    )
    return seen


def test_project_is_resolved_from_ena(monkeypatch):
    seen = _patch(monkeypatch)
    ds = resolve("prjna672658", session=object())

    assert seen == {"runs": "PRJNA672658", "study": "PRJNA672658"}
    assert ds.metadata.archive == "sra"
    assert ds.metadata.resolved_by == "native"
    assert ds.metadata.title == "scRNA of normal and RSA sample"
    assert ds.metadata.relations == {"BioProject": "PRJNA672658", "SRA": "SRP289408"}
    assert ds.metadata.organisms == ["Homo sapiens"]
    assert ds.metadata.contributors == ["Xiamen University"]
    assert len(ds.runs) == 3 and len(ds.files) == 6
    assert {f.md5 for f in ds.files} == {"a" * 32, "b" * 32}


def test_runs_are_grouped_by_biosample(monkeypatch):
    _patch(monkeypatch)
    ds = resolve("PRJNA672658", session=object())

    assert [s.accession for s in ds.samples] == ["SAMN16565585", "SAMN16565584"]
    first = ds.samples[0]
    assert first.source_name == "RSA_5"
    assert first.title == "RSA_5 title"
    assert first.relations["runs"] == "SRR12921178;SRR12921180"
    assert first.characteristics == {"tissue": "decidua", "acc": "SAMN16565585"}
    assert len(first.files) == 4
    assert {r.sample for r in ds.runs} == {"SAMN16565585", "SAMN16565584"}


def test_by_sample_layout_uses_biosample_directories(monkeypatch):
    _patch(monkeypatch)
    ds = resolve("PRJNA672658", session=object())
    apply_layout(ds.files, "by-sample")
    assert "SAMN16565585/SRR12921178_1.fastq.gz" in {f.target for f in ds.files}


def test_study_lookup_uses_the_project_even_for_a_run_accession(monkeypatch):
    seen = _patch(monkeypatch)
    resolve("SRR12921178", session=object())
    assert seen == {"runs": "SRR12921178", "study": "PRJNA672658"}


def test_no_public_runs_gives_an_empty_dataset_not_an_error(monkeypatch):
    _patch(monkeypatch, runs_text=HEADER)
    ds = resolve("PRJNA690251", session=object())
    assert ds.runs == [] and ds.files == [] and ds.samples == []
    assert ds.metadata.title == "scRNA of normal and RSA sample"


def test_include_samples_false_skips_biosamples(monkeypatch):
    _patch(monkeypatch)
    monkeypatch.setattr(
        project, "fetch_characteristics", lambda *a: (_ for _ in ()).throw(AssertionError)
    )
    ds = resolve("PRJNA672658", session=object(), include_samples=False)
    assert ds.samples[0].characteristics == {}


def test_geo_only_options_are_ignored(monkeypatch):
    _patch(monkeypatch)
    ds = resolve(
        "PRJNA672658",
        session=object(),
        include_matrix=True,
        include_archives=True,
        fetch_sizes=False,
        include_reads=False,
    )
    assert len(ds.files) == 6


def test_resolver_choice_is_native_for_projects():
    assert choose_resolver("PRJNA672658", "auto") == "native"
    assert choose_resolver("PRJNA672658", "geofetch") == "native"


def test_raw_ena_records_are_kept(monkeypatch):
    _patch(monkeypatch)
    ds = resolve("PRJNA672658", session=object())
    assert set(ds.raw_documents) == {"PRJNA672658.ena_study.json", "PRJNA672658.ena_runs.tsv"}
    assert ds.raw_documents["PRJNA672658.ena_runs.tsv"].count("\n") == 4


def test_tables(monkeypatch, tmp_path):
    _patch(monkeypatch)
    ds = resolve("PRJNA672658", session=object())
    runs = write_runs_table(ds, tmp_path / "runs.tsv").read_text().splitlines()
    assert runs[0].split("\t")[1] == "sample"
    assert runs[1].split("\t")[:2] == ["SRR12921178", "SAMN16565585"]
    samples = write_samples_table(ds, tmp_path / "samples.tsv").read_text()
    assert "RSA_5" in samples and "decidua" in samples


class _Resp:
    def __init__(self, payload):
        self.payload = payload

    def raise_for_status(self):
        pass

    def json(self):
        return self.payload


class _Session:
    def __init__(self, payload):
        self.payload = payload

    def get(self, url, **kwargs):
        return _Resp(self.payload)


def test_characteristics_drop_insdc_bookkeeping():
    payload = {
        "characteristics": {
            "tissue": [{"text": "Endometrium"}],
            "sex": [{"text": "female"}],
            "INSDC status": [{"text": "live"}],
            "SRA accession": [{"text": "SRS1"}],
            "cell type": [{"text": "T"}, {"text": "NK"}],
        }
    }
    chars = project.fetch_characteristics(_Session(payload), "SAMN1")
    assert chars == {"tissue": "Endometrium", "sex": "female", "cell type": "T; NK"}


def test_characteristics_failure_is_not_fatal():
    class Broken:
        def get(self, url, **kwargs):
            import requests

            raise requests.ConnectionError("down")

    assert project.fetch_characteristics(Broken(), "SAMN1") == {}
