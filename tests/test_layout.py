"""Output layout planning."""

import pytest

from omixfetch.layout import apply_layout, canonical_tenx_name
from omixfetch.models import RemoteFile


def _f(name, accession=None):
    return RemoteFile(url=f"https://x/{name}", filename=name, accession=accession)


@pytest.mark.parametrize(
    "name,expected",
    [
        ("GSM5008737_RNA_3P-barcodes.tsv.gz", "barcodes.tsv.gz"),
        ("GSM5008737_RNA_3P-features.tsv.gz", "features.tsv.gz"),
        ("GSM5008737_RNA_3P-matrix.mtx.gz", "matrix.mtx.gz"),
        ("GSM123_genes.tsv.gz", "genes.tsv.gz"),
        # Not a triplet member.
        ("GSE164378_RAW.tar", None),
        ("GSM123_counts.csv.gz", None),
        ("barcodes.tsv.gz", None),
    ],
)
def test_canonical_tenx_name(name, expected):
    assert canonical_tenx_name(name) == expected


def test_flat_layout_keeps_original_names():
    files = [_f("GSM1_RNA-barcodes.tsv.gz", "GSM1"), _f("GSE1_RAW.tar", "GSE1")]
    apply_layout(files, "flat")
    assert [f.target for f in files] == ["GSM1_RNA-barcodes.tsv.gz", "GSE1_RAW.tar"]


def test_by_sample_groups_and_renames_a_triplet():
    files = [
        _f("GSM5008737_RNA_3P-barcodes.tsv.gz", "GSM5008737"),
        _f("GSM5008737_RNA_3P-features.tsv.gz", "GSM5008737"),
        _f("GSM5008737_RNA_3P-matrix.mtx.gz", "GSM5008737"),
    ]
    apply_layout(files, "by-sample")
    assert [f.target for f in files] == [
        "GSM5008737/barcodes.tsv.gz",
        "GSM5008737/features.tsv.gz",
        "GSM5008737/matrix.mtx.gz",
    ]


def test_by_sample_leaves_series_files_at_top_level():
    files = [_f("GSE1_counts.csv.gz", "GSE1"), _f("GSM1_RNA-matrix.mtx.gz", "GSM1")]
    apply_layout(files, "by-sample")
    assert files[0].target == "GSE1_counts.csv.gz"
    assert files[1].target == "GSM1/matrix.mtx.gz"


def test_collision_within_a_sample_keeps_original_names():
    """A GSM holding two modalities would otherwise overwrite itself."""
    files = [
        _f("GSM1_RNA-matrix.mtx.gz", "GSM1"),
        _f("GSM1_ADT-matrix.mtx.gz", "GSM1"),
    ]
    apply_layout(files, "by-sample")
    assert [f.target for f in files] == [
        "GSM1/GSM1_RNA-matrix.mtx.gz",
        "GSM1/GSM1_ADT-matrix.mtx.gz",
    ]
    # Crucially, the two files still have distinct destinations.
    assert len({f.target for f in files}) == 2


def test_by_sample_does_not_rename_non_triplet_files():
    files = [_f("GSM1_counts.csv.gz", "GSM1")]
    apply_layout(files, "by-sample")
    assert files[0].target == "GSM1/GSM1_counts.csv.gz"


def test_raw_reads_group_under_their_geo_sample():
    """ENA runs carry the GEO sample accession once mapped."""
    files = [_f("SRR1_1.fastq.gz", "GSM1"), _f("SRR1_2.fastq.gz", "GSM1")]
    apply_layout(files, "by-sample")
    assert [f.target for f in files] == ["GSM1/SRR1_1.fastq.gz", "GSM1/SRR1_2.fastq.gz"]


def test_unknown_layout_rejected():
    with pytest.raises(ValueError, match="Unknown layout"):
        apply_layout([_f("a.gz", "GSM1")], "nested")
