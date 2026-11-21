"""GEO filelist.txt parsing and archive deduplication."""

from omixfetch.sources.filelist import archive_members, parse_filelist

FILELIST = (
    "#Archive/File\tName\tTime\tSize\tType\n"
    "Archive\tGSE164378_RAW.tar\t01/17/2021 06:19:10\t1467842560\tTAR\n"
    "File\tGSM5008737_RNA_3P-barcodes.tsv.gz\t01/04/2021 09:05:23\t683095\tTSV\n"
    "File\tGSM5008737_RNA_3P-matrix.mtx.gz\t01/04/2021 09:05:23\t1125353463\tMTX\n"
)


def test_parse_filelist_reads_sizes_and_marks_the_archive():
    entries = parse_filelist(FILELIST)
    assert len(entries) == 3

    archive = entries[0]
    assert archive.name == "GSE164378_RAW.tar"
    assert archive.size == 1467842560
    assert archive.is_archive is True

    member = entries[1]
    assert member.size == 683095
    assert member.type == "TSV"
    assert member.is_archive is False


def test_parse_filelist_tolerates_a_missing_size():
    text = "#Archive/File\tName\tTime\tSize\tType\nFile\tx.gz\t-\t\tTSV\n"
    assert parse_filelist(text)[0].size is None


def test_parse_filelist_empty():
    assert parse_filelist("") == []


def test_archive_members_maps_tar_to_its_contents():
    members = archive_members(parse_filelist(FILELIST))
    assert set(members) == {"GSE164378_RAW.tar"}
    assert members["GSE164378_RAW.tar"] == {
        "GSM5008737_RNA_3P-barcodes.tsv.gz",
        "GSM5008737_RNA_3P-matrix.mtx.gz",
    }


def test_archive_members_empty_without_an_archive_row():
    text = "#Archive/File\tName\tTime\tSize\tType\nFile\tx.gz\t-\t10\tTSV\n"
    assert archive_members(parse_filelist(text)) == {}
