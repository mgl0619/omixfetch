"""Shared fixtures. Every test runs offline against recorded archive records."""

from __future__ import annotations

import pytest

SERIES_SOFT = """^SERIES = GSE000001
!Series_title = A small test series
!Series_geo_accession = GSE000001
!Series_status = Public on Jan 15 2021
!Series_submission_date = Jan 07 2021
!Series_last_update_date = Jan 15 2021
!Series_summary = First summary paragraph.
!Series_summary = Second summary paragraph.
!Series_overall_design = Two samples, one condition.
!Series_type = Expression profiling by high throughput sequencing
!Series_contributor = Ada,,Lovelace
!Series_contributor = Alan,,Turing
!Series_sample_id = GSM000001
!Series_sample_id = GSM000002
!Series_sample_organism = Homo sapiens
!Series_platform_id = GPL24676
!Series_pubmed_id = 33456401
!Series_relation = BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA000001
!Series_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRP000001
!Series_contact_name = Ada,,Lovelace
!Series_contact_institute = Test Institute
!Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE000001/suppl/GSE000001_RAW.tar
!Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE000001/suppl/GSE000001_counts.csv.gz
"""

SAMPLES_SOFT = """^SAMPLE = GSM000001
!Sample_title = Control replicate 1
!Sample_geo_accession = GSM000001
!Sample_source_name_ch1 = PBMC
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = cell type: PBMC
!Sample_characteristics_ch1 = treatment: control
!Sample_library_strategy = RNA-Seq
!Sample_library_source = transcriptomic
!Sample_library_selection = cDNA
!Sample_instrument_model = Illumina NovaSeq 6000
!Sample_relation = BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN00000001
!Sample_supplementary_file_1 = ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnn/GSM000001/suppl/GSM000001_matrix.mtx.gz
^SAMPLE = GSM000002
!Sample_title = Treated replicate 1
!Sample_geo_accession = GSM000002
!Sample_source_name_ch1 = PBMC
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = cell type: PBMC
!Sample_characteristics_ch1 = treatment: stimulated
!Sample_characteristics_ch1 = timepoint: 24h
!Sample_library_strategy = RNA-Seq
!Sample_supplementary_file_1 = ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnn/GSM000002/suppl/GSM000002_matrix.mtx.gz
!Sample_supplementary_file_2 = NONE
"""

@pytest.fixture
def series_soft() -> str:
    return SERIES_SOFT


@pytest.fixture
def samples_soft() -> str:
    return SAMPLES_SOFT
BIOSTUDIES_JSON = {
    "accno": "E-MTAB-0001",
    "attributes": [{"name": "ReleaseDate", "value": "2016-09-21"}],
    "type": "submission",
    "section": {
        "accno": "E-MTAB-0001",
        "type": "Study",
        "attributes": [
            {"name": "Title", "value": "A small test experiment"},
            {"name": "Description", "value": "Single-cell RNA-seq of test tissue."},
            {"name": "Study Type", "value": "RNA-seq of coding RNA from single cells"},
            {"name": "Organism", "value": "Homo sapiens"},
        ],
        "subsections": [
            {
                "type": "Processed Data",
                "files": [
                    [
                        {
                            "path": "counts.txt",
                            "size": 722957491,
                            "type": "file",
                            "attributes": [
                                {"name": "Description", "value": "Processed Data"}
                            ],
                        }
                    ]
                ],
            },
            {
                "type": "MAGE-TAB",
                "files": [
                    {
                        "path": "E-MTAB-0001.idf.txt",
                        "size": 4201,
                        "type": "file",
                        "attributes": [{"name": "Type", "value": "IDF File"}],
                    },
                    {
                        "path": "E-MTAB-0001.sdrf.txt",
                        "size": 2665069,
                        "type": "file",
                        "attributes": [{"name": "Type", "value": "SDRF File"}],
                    },
                ],
            },
        ],
    },
}

IDF_TEXT = """Comment[ArrayExpressAccession]\tE-MTAB-0001
Investigation Title\tA small test experiment
Experiment Description\tSingle-cell RNA-seq of test tissue.
Protocol Name\tP-MTAB-1\tP-MTAB-2
Protocol Description\tCells were dissociated.\tLibraries were prepared.
Public Release Date\t2016-09-21
"""

SDRF_TEXT = (
    "Source Name\tCharacteristics [organism]\tCharacteristics [cell type]\t"
    "Comment[ENA_SAMPLE]\tComment[LIBRARY_STRATEGY]\tAssay Name\n"
    "cell_1\tHomo sapiens\tbeta cell\tERS1234567\tRNA-Seq\tassay_1\n"
    "cell_2\tHomo sapiens\talpha cell\tERS1234568\tRNA-Seq\tassay_2\n"
)




@pytest.fixture
def biostudies_json() -> dict:
    return BIOSTUDIES_JSON


@pytest.fixture
def idf_text() -> str:
    return IDF_TEXT


@pytest.fixture
def sdrf_text() -> str:
    return SDRF_TEXT
