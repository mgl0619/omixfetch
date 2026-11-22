# Roadmap

Planned archives and the groundwork already established for each. Nothing
here is implemented yet beyond what the status column says.

| Archive | Accessions | Status |
|---|---|---|
| NCBI GEO | GSE, GSM, GPL | shipped (0.1–0.5) |
| ENA reads | via GEO BioProject/SRP | shipped (0.2) |
| ArrayExpress / BioStudies | `E-MTAB-*`, `E-GEOD-*` | parsers written, branch `arrayexpress-support` |
| PRIDE (proteomics) | `PXD-*` | API verified, not started |
| MetaboLights (metabolomics) | `MTBLS*` | API partly verified, not started |

The resolver seam makes these additive: a resolver's only contract is to
return a `Dataset` (study metadata, samples, `RemoteFile` list). Everything
downstream — filtering, layout, transfer, manifests, PEP export — is already
archive-agnostic.

---

## ArrayExpress / BioStudies — `E-MTAB-*`

**Status:** resolver and tests written, parked on `arrayexpress-support`.
Needs registering in `resolve.py`, re-enabling `E-*` in `parse_accession`,
plus CLI and docs.

Verified against the live API:

- Study metadata: `https://www.ebi.ac.uk/biostudies/api/v1/studies/{ACC}`.
  Files are nested inside `section.subsections`, sometimes double-nested in
  lists, so collection has to recurse.
- Per-sample annotation comes from the MAGE-TAB SDRF; E-MTAB-5061 has
  3,514 rows × 49 columns. Study-level annotation comes from the IDF.
- Downloads: `https://www.ebi.ac.uk/biostudies/files/{ACC}/{path}` —
  HTTP 200, `Accept-Ranges: bytes`, so aria2 splitting works.
- **No checksums.** File records carry name, size, type, format and
  description only. Verification would be size-only, as with GEO.
- Many E-MTAB experiments link to ENA runs, so the existing `sources/ena.py`
  should serve raw reads once study→project mapping is added.

## PRIDE — `PXD-*`

**Status:** API verified, nothing written.

- File listing: `https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{ACC}/files`
  returns `fileName`, `fileSizeBytes`, `checksum`, `fileCategory` and
  `publicFileLocations`.
- Download: the FTP location rewrites cleanly to HTTPS
  (`https://ftp.pride.ebi.ac.uk/pride/data/archive/...`) and answers
  HTTP 200 with `Accept-Ranges: bytes`. Confirmed on a 1.28 GB `.raw` file.
- An Aspera (`fasp://`) location is also offered and is ignored.
- **Checksums are SHA-1, not md5** (40 hex characters). `RemoteFile.md5`
  would need to become a generic `(algorithm, value)` pair; aria2 accepts
  `checksum=sha-1=...` in its input file, and the fallback backend would
  need the same generalisation.
- Data model differs from GEO: a PXD project has files and assays but no
  "sample → supplementary file" relationship, so the by-sample layout should
  group by `fileCategory` (RAW / PEAK / SEARCH / OTHER) instead.

## MetaboLights — `MTBLS*`

**Status:** listing endpoint verified, download base not yet established.

- File listing: `https://www.ebi.ac.uk/metabolights/ws/studies/{ACC}/files`
  returns entries with `file`, `type` (`metadata_assay`, `metadata_sample`,
  `raw`, ...), `directory` and `status`.
- Study metadata: `https://www.ebi.ac.uk/metabolights/ws/studies/{ACC}`
  returns an ISA-Tab-derived payload.
- **Open question:** the listing gives no size, no checksum and no URL, and
  the obvious FTP path
  (`ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/{ACC}/{file}`)
  returned 404 for MTBLS1. The real download base has to be established
  before anything is built; sizes would otherwise need a HEAD per file.
- Sample annotation lives in ISA-Tab `s_*.txt` and assay `a_*.txt` tables,
  which need a parser similar to the SDRF one already written for MAGE-TAB.

---

---

## Decision: conversion moves to a separate tool

omixfetch's job is **downloading**: finding the right files for an accession
(processed data, raw reads as published, metadata) and fetching them fast and
verified. Converting between formats (`.sra` → FASTQ, BAM → FASTQ,
`.sra` → BAM) belongs in a separate tool with its own dependencies.

Why: conversion drags in external binaries (SRA Toolkit, `bamtofastq`,
`samtools`), platform quirks (the osx-arm64 sra-tools/aria2 conda conflict,
Rosetta for `bamtofastq`, NCBI's reference service for aligned `.sra`) and a
different failure surface from downloading. Keeping them apart lets each be
installed, versioned and tested on its own.

Until the split, the existing conversion code stays as it is, so nothing that
works breaks: `--convert fastq`, `omixfetch fastq` and `omixfetch setup`.
New work goes into downloading. When the split happens:

- **omixfetch keeps:** resolution (GEO, ENA, NCBI SDL), the route knowledge
  of *which files* to fetch, including the BAM-header probe that tells a
  user whether a BAM can later be converted back, plus transfer backends,
  layout, metadata, manifests and the HTML report.
- **The conversion tool takes:** `convert.py`, `bundled.py`, the `fastq`
  pipeline's convert steps, and the conversion half of the decision tree.
  It reads omixfetch's manifest as input, so the two compose.

## Deferred: raw-read conversion

The raw-read side (`omixfetch fastq`) works for the cases tested and makes
the safe choice by default, but these items were put on hold as lower
priority than processed data. Recorded so the findings are not lost.

### Needs a real run on an unsandboxed machine

The development sandbox blocks NCBI's reference service, which every
aligned-`.sra` operation contacts first (`fasterq-dump`, `sam-dump`,
`align-info`, even `vdb-dump` on a local file all hung waiting for it).

| Item | How to test | Status |
|---|---|---|
| 10x BAM → `bamtofastq` | `OMIXFETCH_LIVE=medium python -m pytest tests/live` (GSE220167 / SRR22543607, 8.3 GB) | not run |
| Aligned `.sra` → FASTQ | `fasterq-dump --split-3 SRR390728.sra` (195 MB) | failed in sandbox only |
| Aligned `.sra` → BAM | `sam-dump -r -u SRR390728.sra \| samtools view -b -o out.bam -` | hung in sandbox |
| Non-10x BAM → `samtools fastq` | any bulk BAM submission whose aligner kept unmapped reads | not run |
| Aspera (`--backend aspera`) | needs `ascp` and outbound UDP 33001 | not run |
| `s5cmd` backend | `--backend s5cmd` against an NCBI S3 route | unit-tested only |

### Possible additions

- **`.sra` → BAM route** via `sam-dump`, for bulk aligned submissions where
  NCBI kept no original BAM (older runs such as SRR390728). Not useful for
  single-cell: a BAM rebuilt from `.sra` lacks the submitter's tags, so a 10x
  run would still have no barcodes. Needs reference access at conversion.
- **SRA Lite** (`.lite`), which NCBI serves at about half the size of the
  normalised `.sra` (22.2 GB against 48.1 GB for SRR16036138) with quality
  scores reduced to pass/fail. Reasonable as an opt-in for bulk
  quantification. Wrong as a default for single-cell, where Cell Ranger and
  STARsolo use quality scores for barcode correction. Served from NCBI rather
  than S3.
- **Bundling `samtools`**, which has no official prebuilt binaries;
  currently it must be installed separately.

### Findings to keep

- NCBI serves three file types per run: the normalised `.sra` (S3), SRA Lite
  (NCBI servers), and the original submission (S3, only for non-FASTQ
  submissions). Confirmed with SDL `filetype=sra|sralite`.
- For a 10x BAM submission, the `.sra`, SRA Lite and ENA's FASTQ all lack the
  barcode read; only the original BAM can be converted back to complete
  FASTQ. GSE184708 is the documented example (see README).
- Before fetching raw data, check GEO's processed files. GSE184708 ships a
  1.2 GiB count matrix against 1.7 TiB of BAMs.

---

## Cross-cutting work these imply

1. **Generalise checksums.** `RemoteFile.md5` assumes one algorithm. PRIDE
   publishes SHA-1, ENA md5, GEO and BioStudies none. Replace with an
   algorithm-tagged field and teach both backends to use it.
2. **Generalise the layout.** `layout.py` hardcodes `GSM\d+` for sample
   grouping and 10x triplet naming. Each archive needs its own notion of
   "group files by", with the GEO behaviour as one case.
3. **Accession routing.** `parse_accession` currently rejects everything
   non-GEO with a GEO-specific message; it needs a registry keyed by pattern.
4. **Decide on metadata harmonisation.** `SampleMetadata` is GEO-shaped.
   SDRF, ISA-Tab and PRIDE assay metadata do not map onto it cleanly, and
   forcing them to may lose more than it gains — the raw documents are
   already preserved alongside the parsed form.
