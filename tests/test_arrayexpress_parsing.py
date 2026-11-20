"""BioStudies JSON, IDF and SDRF parsing."""

from omixfetch.sources.arrayexpress import (
    _attributes,
    _sdrf_to_samples,
    _walk_files,
    file_url,
    parse_idf,
    parse_sdrf,
)


def test_walk_files_finds_nested_and_doubly_nested_entries(biostudies_json):
    collected = []
    _walk_files(biostudies_json["section"], collected)
    paths = {c["path"] for c in collected}
    # counts.txt is wrapped in a nested list; the MAGE-TAB files are not.
    assert paths == {"counts.txt", "E-MTAB-0001.idf.txt", "E-MTAB-0001.sdrf.txt"}


def test_attributes_joins_repeated_names():
    node = {"attributes": [{"name": "X", "value": "1"}, {"name": "X", "value": "2"}]}
    assert _attributes(node) == {"X": "1; 2"}


def test_file_url():
    assert file_url("E-MTAB-0001", "counts.txt") == (
        "https://www.ebi.ac.uk/biostudies/files/E-MTAB-0001/counts.txt"
    )


def test_parse_idf_collects_multivalued_rows(idf_text):
    idf = parse_idf(idf_text)
    assert idf["Investigation Title"] == ["A small test experiment"]
    assert idf["Protocol Name"] == ["P-MTAB-1", "P-MTAB-2"]
    assert len(idf["Protocol Description"]) == 2


def test_parse_sdrf_rows_and_duplicate_headers(sdrf_text):
    rows = parse_sdrf(sdrf_text)
    assert len(rows) == 2
    assert rows[0]["Source Name"] == "cell_1"
    assert rows[1]["Characteristics [cell type]"] == "alpha cell"


def test_parse_sdrf_suffixes_repeated_columns():
    text = "A\tA\tB\n1\t2\t3\n"
    rows = parse_sdrf(text)
    assert rows[0] == {"A": "1", "A#2": "2", "B": "3"}


def test_sdrf_to_samples_extracts_characteristics_and_relations(sdrf_text):
    samples = _sdrf_to_samples(parse_sdrf(sdrf_text))
    assert len(samples) == 2
    first = samples[0]
    assert first.accession == "cell_1"
    assert first.organism == "Homo sapiens"
    assert first.characteristics == {"organism": "Homo sapiens", "cell type": "beta cell"}
    assert first.relations["ENA_SAMPLE"] == "ERS1234567"
    assert first.library_strategy == "RNA-Seq"
