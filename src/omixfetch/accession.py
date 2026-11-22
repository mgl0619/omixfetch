"""Accession recognition and routing."""

from __future__ import annotations

import enum
import re

from .errors import AccessionError


class Archive(str, enum.Enum):
    """Archives this release can resolve."""

    GEO = "geo"
    #: Sequencing projects in the INSDC read archives (SRA, ENA, DDBJ),
    #: addressed by BioProject, study, experiment, run or BioSample accession.
    SRA = "sra"


#: GEO series, sample and platform accessions.
GEO_RE = re.compile(r"^(GSE|GSM|GPL)(\d+)$", re.IGNORECASE)

#: Read-archive accessions ENA's portal resolves to runs: BioProjects
#: (PRJNA/PRJEB/PRJDB), studies (SRP/ERP/DRP), experiments (SRX/ERX/DRX),
#: runs (SRR/ERR/DRR) and BioSamples (SAMN/SAMEA/SAMD). Whichever archive
#: received the submission, all three mirror it, so one pattern covers them.
SRA_RE = re.compile(
    r"^(PRJ(?:NA|EB|DB)|[SED]RP|[SED]RX|[SED]RR|SAMN|SAMEA|SAMD)(\d+)$",
    re.IGNORECASE,
)

#: ArrayExpress accessions. Recognised so the error message can be specific,
#: but not resolvable until ArrayExpress support lands (see ROADMAP).
ARRAYEXPRESS_RE = re.compile(r"^E-[A-Z]{4}-\d+$", re.IGNORECASE)


def parse_accession(accession: str) -> tuple[Archive, str]:
    """Classify *accession* and return ``(archive, normalised_accession)``.

    Accessions are normalised to canonical upper case, so ``gse164378`` and
    ``GSE164378`` resolve identically.

    Raises:
        AccessionError: if the accession is malformed, or belongs to an
            archive this release does not yet support.
    """
    acc = accession.strip()
    if not acc:
        raise AccessionError("Empty accession.")

    if GEO_RE.match(acc):
        return Archive.GEO, acc.upper()

    if SRA_RE.match(acc):
        return Archive.SRA, acc.upper()

    if ARRAYEXPRESS_RE.match(acc):
        raise AccessionError(
            f"{acc.upper()} is an ArrayExpress accession. This release of omixfetch "
            f"covers NCBI GEO and the SRA/ENA read archives; ArrayExpress support is in progress "
            f"on the 'arrayexpress-support' branch and is planned for 0.2."
        )

    raise AccessionError(
        f"Unrecognised accession {accession!r}. Supported: GEO series, samples and "
        f"platforms (GSE/GSM/GPL), and read-archive projects, studies, experiments, "
        f"runs and samples (PRJNA/PRJEB/PRJDB, SRP/ERP/DRP, SRX/ERX/DRX, "
        f"SRR/ERR/DRR, SAMN/SAMEA/SAMD)."
    )


def geo_ftp_stub(accession: str) -> str:
    """Return the GEO FTP directory stub for *accession*.

    GEO buckets records in directories whose final three digits are replaced
    by ``nnn`` -- ``GSE164378`` lives under ``GSE164nnn``. Accessions with
    three or fewer digits live under a bare ``nnn`` bucket.
    """
    m = GEO_RE.match(accession)
    if not m:
        raise AccessionError(f"{accession!r} is not a GEO accession.")
    prefix, digits = m.group(1).upper(), m.group(2)
    bucket = digits[:-3] + "nnn" if len(digits) > 3 else "nnn"
    return f"{prefix}{bucket}"


def geo_kind(accession: str) -> str:
    """Return the GEO FTP sub-tree (``series``/``samples``/``platforms``)."""
    m = GEO_RE.match(accession)
    if not m:
        raise AccessionError(f"{accession!r} is not a GEO accession.")
    return {"GSE": "series", "GSM": "samples", "GPL": "platforms"}[m.group(1).upper()]
