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
