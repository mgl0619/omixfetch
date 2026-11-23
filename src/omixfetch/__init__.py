"""omixfetch: parallel downloads of public genomics datasets.

This release resolves NCBI GEO accessions into a file manifest plus
harmonised study and sample metadata, optionally including raw sequencing
reads mirrored by ENA, then hands the transfer to aria2c for
multi-connection, resumable downloading.
"""

__version__ = "0.9.0"

from .accession import Archive, parse_accession
from .models import Dataset, DatasetMetadata, RemoteFile, SampleMetadata
from .resolve import resolve

__all__ = [
    "__version__",
    "RemoteFile",
    "SampleMetadata",
    "DatasetMetadata",
    "Dataset",
    "parse_accession",
    "Archive",
    "resolve",
]
