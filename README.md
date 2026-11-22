# omixfetch

Parallel, resumable downloads of public genomics datasets from **NCBI GEO** --
processed files and raw sequencing reads -- with full study and sample metadata
extracted alongside the data.

`omixfetch` resolves an accession into a file manifest plus harmonised metadata,
then hands the transfer to [aria2c](https://aria2.github.io/) — which already
solves multi-connection range requests, correct resumption of partial
downloads, retry with backoff, and checksum verification. The Python side does
the part no existing tool does well: turning an accession into a complete,
filterable manifest and a tidy metadata table.

## Install

**Python 3.10 or newer is required.** Python 3.9 is past end of life and is not
supported; `pip install` will refuse it rather than fail at import time.

### With mamba or conda (recommended)

`environment.yml` pins a working interpreter and includes `aria2`, so one
command gives a complete setup:

```bash
mamba env create -f environment.yml     # or: micromamba / conda
mamba activate omixfetch
pip install -e ".[dev]"
pytest
```

Using `micromamba` for the first time? It needs a root prefix before
`env create` works:

```bash
micromamba shell init -s zsh -p ~/micromamba   # then restart your shell
```

### With a plain virtualenv

No conda needed, but install `aria2` separately:

```bash
python3.13 -m venv .venv                 # any interpreter >= 3.10
source .venv/bin/activate
pip install -e ".[dev]"
brew install aria2                       # macOS; apt install aria2 on Debian/Ubuntu
pytest
```

### From PyPI

```bash
pip install omixfetch
```

`aria2c` is strongly recommended and used automatically when present:

```bash
brew install aria2                        # macOS
sudo apt install aria2                    # Debian/Ubuntu
conda install -c conda-forge aria2        # conda
```

> **If `pytest` reports `ModuleNotFoundError: No module named 'omixfetch'`**, the
> editable install did not take effect — usually because `pip` and `pytest` are
> different interpreters, or because the active Python is 3.9. Check with
> `python -c "import sys; print(sys.version)"` and `python -m pytest` (which
> always uses the interpreter you think it does).

Without it, omixfetch falls back to a threaded pure-Python downloader that
transfers several files at once but does not split a single file across
connections. The same fallback kicks in if `aria2c` is installed but cannot
complete a transfer -- some environments put an intercepting TLS proxy in
front of the network that aria2's TLS backend refuses. Force either engine
with `--backend aria2` or `--backend python`.

## Usage

Inspect before committing to a download:

```bash
omixfetch info GSE164378
omixfetch list GSE164378 --include "*.csv.gz"
```

Download data and metadata, 4 files at a time with 16 connections each:

```bash
omixfetch get GSE164378 -o ./data -j 4 -x 16
```

Useful options:

```bash
omixfetch get GSE164378 -o ./data \
    --include "*.h5" --include "*.mtx.gz" \
    --exclude "*_RAW.tar" \
    --max-size 5GB \
    --matrix            # also fetch the series matrix file
omixfetch get GSE164378 --dry-run          # resolve and report, download nothing
omixfetch metadata GSE164378 -o ./data     # metadata only, no data files
```

### Metadata output and PEP

Every download writes metadata beside the data. Add `--pep` to additionally
export a [PEP](https://pep.databio.org/) project — a `config.yaml` plus a
sample CSV — which `peppy`, `looper` and `eido` read directly:

```bash
omixfetch metadata GSE164378 -o ./data --pep
```

```yaml
pep_version: 2.1.0
name: GSE164378
sample_table: GSE164378_samples.csv
experiment_metadata:
  title: Integrated analysis of multimodal single-cell data
  organisms: [Homo sapiens]
```

PEP is written natively, so this costs one small dependency (`pyyaml`) rather
than the full pepkit stack.

### Choosing a resolver

Metadata resolution has two backends, selected by `--resolver`:

| Value | Behaviour |
|---|---|
| `auto` (default) | geofetch for GSE series, except where it cannot serve the request |
| `native` | the built-in SOFT parser, always |
| `geofetch` | geofetch, always (GSE only) |

`auto` prefers [geofetch](https://pep.databio.org/geofetch/) because it is the
more mature and widely exercised GEO parser:

```bash
pip install "omixfetch[geofetch]"     # auto then uses it for GSE
omixfetch info GSE164378              # the Resolver row shows which one ran
```

`auto` falls back to the native resolver in three cases, each verified against
geofetch 0.12.11:

- **GSM and GPL accessions** — geofetch's project API returns an empty result
  for GSM and raises on GPL.
- **Raw-read workflows** (`--reads`, `--reads-only`) — geofetch's sample table
  has no `relation` columns, so the BioSample identifiers linking an ENA run to
  its GEO sample are unavailable and runs would be left ungrouped. The native
  resolver maps all 16 runs of GSE52778 to their samples.
- **Series with no processed files** — geofetch returns a series table and no
  samples at all, which is common for reads-only submissions.

For processed-file workflows both agree on GSE164378: 27 files, 1.47 GB.

**Transfers go through aria2c either way.** geofetch is called with
`just_metadata=True` and never downloads anything; it resolves URLs, and
omixfetch moves the bytes.

**The SRA Toolkit is not required by either path.** `pip install geofetch`
installs Python packages only — `prefetch` and `fasterq-dump` are separate
binaries that geofetch shells out to for its own raw-data path. omixfetch
fetches reads from ENA instead, so you never need them. If you do want them:
`mamba install -c bioconda sra-tools`.

### Raw reads

`--reads` additionally resolves FASTQ files through the ENA portal, which
mirrors SRA and serves reads as ordinary HTTPS objects. No SRA Toolkit, and no
`.sra` conversion:

```bash
omixfetch get GSE52778 -o ./data --reads          # processed files and FASTQ
omixfetch get GSE52778 -o ./data --reads-only     # FASTQ only
omixfetch list GSE52778 --reads                   # 39 files, 46.9 GB
```

ENA publishes a **per-file md5**, so raw reads are verified end to end — both
backends check the checksum and refuse a corrupt file rather than leaving it on
disk. GEO publishes no checksums for supplementary files, so those are verified
on declared size only.

#### Choosing where the bytes come from

The same reads are served by EBI and by NCBI, and which is faster depends
entirely on your network path. `--source` selects:

```bash
omixfetch get GSE184708 -o ./data --reads-only --source ena    # EBI (default)
omixfetch get GSE184708 -o ./data --reads-only --source ncbi   # NCBI SDL / S3
```

`--source ncbi` resolves locations through NCBI's SDL service — the same one
the SRA Toolkit queries, but plain HTTP, so no toolkit is needed to download.
It exposes two objects per run:

| `--sra-flavour` | What you get |
|---|---|
| `sra` (default) | the normalised `.sra` archive — smaller, convert with `fasterq-dump` |
| `original` | the file the submitter uploaded, often a BAM |

For GSE184708 the `.sra` set is 709 GiB against 1.7 TiB of submitted BAMs.

#### Download and convert in one step

If you want FASTQ, ask for it directly:

```bash
omixfetch get GSE184708 -o ./data --reads-only --source ncbi --convert fastq -t 8
```

This downloads each `.sra`, runs `fasterq-dump`, gzips the result and removes
the archive (`--keep-sra` to retain it).

Conversion needs the SRA Toolkit, which pip does not install. **It does not
have to live in omixfetch's environment** — on osx-arm64 `sra-tools` and
`aria2` cannot be solved into one environment at all: they disagree over
`libxml2`/`icu`, and a combined solve fails outright. So the toolkit gets its
own:

```bash
mamba env create -f environment-convert.yml
export OMIXFETCH_FASTERQ_DUMP="$(conda run -n omixfetch-convert which fasterq-dump)"
```

Or point at any binary directly, including one from a container or NCBI's
dependency-free [release tarballs](https://github.com/ncbi/sra-tools/releases):

```bash
omixfetch get ... --convert fastq --fasterq-dump /opt/sra-tools/bin/fasterq-dump
```

Downloading `.sra` archives needs no toolkit at all; only conversion does.

#### Measuring which source is fastest

Published throughput figures do not transfer between networks. Measure yours:

```bash
omixfetch benchmark GSE184708 --sample 100MB
```

It samples a fixed byte range from each source and reports MiB/s plus the
projected time for the whole dataset, costing a few hundred MB rather than
terabytes.

#### Aspera

ENA publishes an Aspera endpoint for most runs. `ascp` often beats HTTPS
substantially over long paths:

```bash
conda install -c hcc aspera-cli
omixfetch get GSE184708 -o ./data --reads-only --backend aspera
```

Aspera needs outbound **UDP port 33001**, which many institutional firewalls
block; the backend says so rather than hanging.

#### FASTQ, BAM and CRAM

ENA serves whatever the submitter deposited. `--prefer` chooses between
ENA-generated FASTQ and the submitter's own files, which for long-read data are
usually BAM:

```bash
omixfetch get GSE52778 -o ./data --reads                      # FASTQ (default)
omixfetch get GSE52778 -o ./data --reads --prefer submitted   # BAM/CRAM if present
```

Either flavour falls back to the other when absent, so a run with only BAMs
still resolves. Format is irrelevant to the transfer — aria2c splits a 12.6 GB
BAM across connections exactly as it does a 2.4 GB FASTQ, and both carry an md5.

Runs are matched back to their GEO samples through the BioSample accession and
written to `<ACC>.runs.tsv`.

**Controlled access.** Many human series have a BioProject but no public runs,
because the reads are restricted through dbGaP. omixfetch says so explicitly
rather than reporting an empty download:

```
GSE164378: ENA serves no public runs for this project. Human data under
controlled access (dbGaP) is the usual reason ...
```

No tool can retrieve those without authorised access.

Emit a manifest and run the transfer yourself:

```bash
omixfetch manifest GSE164378 -o files.txt --format aria2
aria2c -i files.txt -j 8 -x 16
```

## Scope: downloading

omixfetch's job is downloading: the right files for an accession, fetched
fast and verified. Format conversion (`.sra` or BAM to FASTQ) is currently
included through `--convert fastq`, `omixfetch fastq` and `omixfetch setup`,
but is slated to move into a separate tool; see [ROADMAP.md](ROADMAP.md).
New features go into downloading.

## Download routes

Every way to fetch one GEO series, in one place. The examples use GSE184708, a
20-run mouse scRNA-seq series submitted as aligned BAMs (1.7 TiB). **Treat it
as a final-stage test, not something to iterate on:** a full fetch takes many
hours.

### Hands-off: the autonomous pipeline

```bash
omixfetch fastq GSE184708 -o ./data --dry-run   # show the plan, fetch nothing
omixfetch fastq GSE184708 -o ./data             # fetch every run as FASTQ
```

The converters it needs — NCBI's SRA Toolkit and 10x `bamtofastq` — are
installed automatically on first use, as prebuilt binaries under
`~/.omixfetch/tools` (set `OMIXFETCH_HOME` to put them elsewhere, e.g. on a
cluster with a small home quota). No conda is involved, so the sra-tools/aria2
conflict does not arise. To install them ahead of time: `omixfetch setup`.

#### How it decides: the decision tree

By default (`--routes auto`) each run follows this tree. It's driven by what
the submitter **actually uploaded**, because that fixes what can be recovered:

```
run
├─ no public files ............................ skip (controlled access)
├─ submitted as BAM/CRAM   (header read with a 4 MB range request)
│   ├─ 10x_bam_to_fastq lines ................. AWS BAM → bamtofastq
│   │                                            .sra and ENA FASTQ excluded: derived
│   │                                            from the BAM, they lost the barcode read
│   ├─ single-cell, no 10x lines .............. no faithful route
│   ├─ aligner kept unmapped reads ............ AWS BAM → samtools fastq
│   └─ dropped unmapped reads / unknown ....... no faithful route
└─ submitted as FASTQ
    └─ ENA FASTQ  or  AWS .sra → fasterq-dump
       ranked by measured throughput × size, conversion included
```

The submission type comes from NCBI's record of the original file, backed by
ENA. Single-cell status comes from the BAM header, because NCBI's
`LibrarySource` can't be trusted: it says `TRANSCRIPTOMIC` for GSE184708,
which is 10x data. Throughput is measured once per series (about 16 MB per
source; `--no-measure` to skip), so "fastest" means fastest *from your
machine*. On a network where ENA is slow, a FASTQ submission is fetched as
`.sra` from AWS and converted, rather than as ENA's ready-made FASTQ.

"No faithful route" fails the run with the reason rather than producing
incomplete FASTQ. `--allow-incomplete` accepts the loss. `--routes sra,original`
forces an order, but data-loss exclusions still apply.

For single-cell `.sra` conversions, `fasterq-dump` runs with
`--split-files --include-technical`. SRA often flags barcode and index reads as
technical, and the usual `--skip-technical` would drop them silently.

#### The report

Every run writes `<ACC>/download_report.html`: one self-contained file, no
scripts or external assets. It shows each run's submission type, the branch
taken, which routes were not used and why, the estimate, the outcome, spot
counts expected against written, bytes, time and throughput, plus any
failures, the measured rates and the tool versions. A `--dry-run` writes the
same report as a plan, before anything is downloaded.

Useful options: `--routes sra,original` to change the order, `--parallel-runs`,
`--include 'SRR160361*'` to fetch a subset, `--keep-intermediate` to keep the
`.sra`/BAM files.

### By hand: each route individually

| # | What | From | Converter |
|---|---|---|---|
| 1 | Processed files (matrices, tables) | GEO | — |
| 2 | FASTQ | ENA | — |
| 3 | Original submission (e.g. BAM) | ENA | `bamtofastq` |
| 4 | Original submission (e.g. BAM) | NCBI's AWS mirror | `bamtofastq` |
| 5 | Normalised `.sra` | NCBI's AWS mirror | `fasterq-dump` |

```bash
# 1  processed files
omixfetch get GSE184708 -o ./data

# 2  FASTQ, when ENA has it (GSE184708 does not -- it was submitted as BAM)
omixfetch get GSE52778 -o ./data --reads-only

# 3  original BAMs from ENA
omixfetch get GSE184708 -o ./data --reads-only --prefer submitted
omixfetch get GSE184708 -o ./data --reads-only --prefer submitted --backend aspera

# 4  original BAMs from NCBI's AWS mirror
omixfetch get GSE184708 -o ./data --reads-only --source ncbi --sra-flavour original
omixfetch get GSE184708 -o ./data --reads-only --source ncbi --sra-flavour original --backend s5cmd

# 5  .sra from AWS, converted to FASTQ in the same command
omixfetch get GSE184708 -o ./data --reads-only --source ncbi --convert fastq -t 8

# which source is fastest from this machine
omixfetch benchmark GSE184708 --sample 100MB
```

Add `--dry-run` to any `get` to see what would be fetched, and
`--include 'SRR16036138*'` to try a single run first.

### Transfer engines

| `--backend` | Works with | Install |
|---|---|---|
| `aria2` | any HTTPS route | in `environment.yml` |
| `s5cmd` | NCBI AWS routes (4, 5) | `conda install -c conda-forge s5cmd` |
| `aspera` | ENA routes (2, 3); needs outbound UDP 33001 | `conda install -c hcc aspera-cli` |
| `python` | any HTTPS route; one connection per file | built in |

`auto` (the default) picks `s5cmd` when every file is on S3 and s5cmd is
installed, otherwise `aria2`, falling back to `python` if aria2 fails.

### Which to choose

The archives' sizes are fixed; their speeds are not — they depend on your
network path, so run `benchmark` from the machine that will do the work.
For aligned submissions like GSE184708, route 5 moves the least data
(709 GiB of `.sra` against 1.7 TiB of BAM) and lands directly as FASTQ, but its
conversion fetches reference sequences from NCBI (see TUTORIAL.md,
troubleshooting). If that is blocked, use route 4.

### What has been verified

Everything below was run against the live archives, not mocked. Offline unit
tests (203) cover the logic; the live tiers are described in
[TESTING.md](TESTING.md).

| Piece | Status | Evidence |
|---|---|---|
| GEO processed files | ✓ verified | real downloads; sizes and gzip integrity checked |
| ENA FASTQ route | ✓ verified | GSE204978, 2 runs: file sizes match ENA to the byte, md5 checked |
| `.sra` route + `fasterq-dump` | ✓ verified | GSE204978, 2 runs: 23,959 and 26,554 reads in both mates, matching NCBI's spot counts |
| Decision tree, FASTQ branch | ✓ verified | GSE204978: measured NCBI 11.9 vs ENA 0.5 MiB/s, so it picked `.sra` over ENA's ready-made FASTQ |
| Decision tree, 10x BAM branch | ✓ planned correctly | GSE184708 dry run: all 20 runs routed AWS BAM → `bamtofastq` |
| Resume after interruption | ✓ verified | rerun skipped both finished runs; nothing re-downloaded |
| Controlled-access series | ✓ verified | GSE164378 reported as having no public runs; exit code 1 |
| ENA submitted file | ✓ verified | one `.bai` from GSE184708; md5 matched ENA |
| BAM header probe | ✓ verified | GSE184708 header read with a 4 MB range request; `10x_bam_to_fastq` lines found |
| Bundled converters | ✓ verified | SRA Toolkit 3.4.1 and `bamtofastq` 1.4.1 installed in 6 s on macOS arm64 |
| HTML report | ✓ verified | written for a completed run and for a dry-run plan |
| `bamtofastq` on a real BAM | **not yet run** | medium live tier: GSE220167 / SRR22543607, 8.3 GB |
| `samtools fastq` on a real BAM | **not yet run** | needs `samtools` (`conda install -c bioconda samtools`); not bundled |
| Aspera (`ascp`) | **not yet run** | endpoints exist for all 20 runs of GSE184708 |
| `s5cmd` backend | **unit-tested only** | not yet run against S3 |
| Aligned `.sra` conversion | ✗ failed in test | SRR390728 needs NCBI's reference service at conversion time (see below) |

Throughput figures quoted anywhere in this README were measured from one
machine and **do not transfer to yours**. That is why the pipeline measures
on every run, and why `omixfetch benchmark` exists.

### Why the original BAM, not the `.sra`, for 10x submissions

This is the most consequential decision the tool makes, so the evidence is
recorded here. GSE184708 was submitted as 10x Cell Ranger BAMs.

**The `.sra` lost the barcode read.** NCBI's RunInfo for SRR16036138:

| Field | Value | Meaning |
|---|---|---|
| `LibraryLayout` | PAIRED | the library has two reads per fragment |
| `spots_with_mates` | **0** | none of them have a mate in the `.sra` |
| `avgLength` | **98** | 98 bases per spot: the cDNA read only |

A 10x library is R1 (cell barcode + UMI, 26–28 bp) plus R2 (cDNA). The
barcode and UMI were stored as BAM *tags*, not as a read, and SRA's loader
did not keep them. Converting that `.sra` produces FASTQ with no cell
barcodes, which Cell Ranger, STARsolo and kallisto|bustools cannot use. ENA's
FASTQ for BAM submissions is generated from the same data and has the same
gap.

**The BAM still has everything.** Its header, read with a 4 MB range request:

```
@CO  10x_bam_to_fastq:R1(CR:CY,UR:UY)   barcode + UMI, from tags
@CO  10x_bam_to_fastq:R2(SEQ:QUAL)      cDNA
@CO  10x_bam_to_fastq:I1(BC:QT)         sample index
@PG  ID:STAR  VN:STAR_2.5.1b  CL:... --outSAMunmapped Within ...
```

`CR`/`UR` are the *uncorrected* barcode and UMI, so the rebuilt R1 matches
what was sequenced. `--outSAMunmapped Within` means unmapped reads were kept.
The STAR input path (`.../SC_RNA_COUNTER_CS/...`) identifies Cell Ranger's
own pipeline. So `bamtofastq` can recover the complete original FASTQ, at
the cost of a larger transfer: 1.7 TiB of BAM against 709 GiB of `.sra`.

**NCBI's labels are not enough on their own.** For GSE184708, NCBI's
`LibrarySource` says `TRANSCRIPTOMIC`, not single-cell, and only the original
file is tagged `TenX`. Conversely, GSE217662's originals are tagged `TenX`
but their BAM headers have no `10x_bam_to_fastq` lines, so `bamtofastq`
would fail on them after a 280 GB download. The header probe catches both
cases before anything large is fetched.

**The general rule** used by the decision tree: a BAM converts back to the
FASTQ that was sequenced only if it kept the unmapped reads and, for
single-cell data, the barcode and UMI. Most Cell Ranger BAMs do. STAR with
default settings, or HISAT2/Bowtie2 with `--no-unal`, drop unmapped reads,
and converting them gives back only what mapped.

### Before downloading terabytes, check what GEO already has

GSE184708 also ships its processed data on GEO, 1.2 GiB in total:

| File | Size |
|---|---|
| `GSE184708_raw_counts_gonad_all_samples_XX_XY_E10_to_E16.mtx.gz` | 1.2 GiB |
| `GSE184708_mayere_barcodes.tsv.gz` | 431 KiB |
| `GSE184708_mayere_genes.tsv.gz` | 106 KiB |
| `GSE184708_mayere_meta_data.csv.gz` | 7.0 MiB |

If you'll analyse the authors' counts as they are, `omixfetch get GSE184708`
fetches those in minutes. Raw FASTQ (`omixfetch fastq`) is worth its 1.7 TiB
only for reprocessing: a new reference, a different aligner, ambient-RNA
correction from raw droplets. For tools that read the alignment directly
(velocyto, genotype demultiplexing), the BAMs alone suffice and conversion
can be skipped.

### Known issues and caveats

- **Aligned `.sra` needs NCBI's reference service.** An `.sra` loaded from an
  aligned submission stores reads as differences against a reference genome;
  `fasterq-dump` fetches those references at conversion time. Installing
  sra-tools locally does not change this. SRR390728 failed this way in
  testing. References are cached, so it is a one-time cost per study; point
  the cache at scratch with `vdb-config --interactive`. The decision tree
  avoids this case for 10x data anyway.
- **sra-tools and aria2 cannot share a conda environment on osx-arm64**: a
  combined solve fails on a `libxml2`/`icu` conflict. The bundled prebuilt
  toolkit (`omixfetch setup`) sidesteps conda entirely.
- **`bamtofastq` on Apple Silicon** is an x86_64 binary and needs Rosetta 2.
  It runs natively on Linux.
- **Fixed in 0.8.0:** single-cell `.sra` conversions ran with
  `--skip-technical`, which would silently drop barcode and index reads that
  SRA flags as technical. They now use `--split-files --include-technical`.
- **Fixed in 0.7.0:** expected spot counts came from ENA, which reports 0 for
  submissions it never processed into reads, so every BAM-submitted run would
  have failed verification. Counts now come from NCBI RunInfo.

### Example reports

[A completed two-run download of GSE204978](docs/examples/GSE204978_download_report.html)
and [the 20-run plan for GSE184708](docs/examples/GSE184708_plan_report.html)
show what `download_report.html` looks like in practice (download the file
and open it in a browser; GitHub shows HTML as source): one row
per run with submission type, decision branch, routes not used and why,
estimate, outcome, spot counts expected against written, size, time and
throughput.

## What you get

Files are grouped by sample by default, with 10x triplets given their
Cell Ranger names so `scanpy.read_10x_mtx()` and Seurat's `Read10X()` can read
a directory straight off:

```
data/GSE164378/
├── data/
│   ├── GSM5008738/
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── GSM5008739/
│   │   └── ...
│   └── GSE164378_sc.meta.data_3P.csv.gz     # series-level files stay at top
└── metadata/
    ├── GSE164378.metadata.json      # complete parsed record
    ├── GSE164378.samples.tsv        # one row per sample, characteristics as columns
    ├── GSE164378.files.tsv          # file manifest with sizes and URLs
    ├── GSE164378.runs.tsv           # sequencing runs, with --reads
    └── pep/                         # with --pep
        ├── GSE164378_config.yaml
        └── GSE164378_samples.csv
    ├── GSE164378.soft.txt           # unmodified series SOFT record
    └── GSE164378.samples.soft.txt   # unmodified sample SOFT records
```

Use `--layout flat` to keep everything in one directory under the names GEO
published. Renaming is applied only when unambiguous: if one GSM holds several
modalities whose files would collide, the original names are kept for that
sample rather than risk overwriting data.

`samples.tsv` flattens each sample's characteristics into columns, taking the
union across samples so sparsely annotated samples carry empty cells rather
than being dropped. The raw SOFT records are kept so nothing is lost to the
parser's opinions.

## Python API

```python
from omixfetch import resolve
from omixfetch.backends import Aria2Backend

dataset = resolve("GSE164378")
print(dataset.metadata.title, dataset.metadata.n_samples, dataset.total_bytes)

big = [f for f in dataset.files if (f.size or 0) > 1e9]
Aria2Backend(connections=16, concurrent=4).download(big, "./data")
```

## Scope and roadmap

Planned archives — ArrayExpress/E-MTAB, PRIDE proteomics, MetaboLights
metabolomics — and the API groundwork established for each are tracked in
[ROADMAP.md](ROADMAP.md).

**0.1 — GEO processed files.** Series (GSE), samples (GSM) and platforms (GPL):
supplementary files, series matrix files, study metadata, per-sample
characteristics.

**0.2 — Raw reads.** FASTQ via the ENA portal, with md5 verification and runs
mapped back to GEO samples. Controlled-access projects reported explicitly.

**0.3 — Resolver choice and PEP.** Optional geofetch resolver, native PEP
export, file sizes read from GEO's `filelist.txt` instead of one HEAD request
per file, and redundant `_RAW.tar` archives dropped when their members are
listed separately.

**0.4 — ArrayExpress.** E-* accessions via the EBI BioStudies API, with
per-sample annotation parsed from MAGE-TAB SDRF files. Work in progress on the
`arrayexpress-support` branch.

Reads are fetched from ENA rather than through
[SRA Toolkit](https://github.com/ncbi/sra-tools). Where a submitter deposited
files SRA serves but ENA does not mirror, the toolkit remains the right tool;
omixfetch records the SRA and BioProject identifiers in its metadata so you can
hand them straight over.

## Related tools

[GEOfetch](https://pep.databio.org/geofetch/) (GEO + SRA, PEP metadata output),
[GEOquery](https://bioconductor.org/packages/GEOquery/) (R, the most complete
GEO client), [GEOparse](https://geoparse.readthedocs.io/) (Python, SOFT
parsing), [ffq](https://github.com/pachterlab/ffq) (metadata and links only) and
[fastq-dl](https://github.com/rpetit3/fastq-dl) (ENA/SRA reads).
omixfetch differs in covering processed files, raw reads and metadata behind one
accession-driven interface, and in treating the transfer itself as the problem
worth solving.

## Tutorial

A step-by-step walkthrough, from installation to a loaded AnnData, is in
[TUTORIAL.md](TUTORIAL.md).

## Development

```bash
mamba env create -f environment.yml && mamba activate omixfetch
pip install -e ".[dev]"
python -m pytest
ruff check src tests
```

The test suite runs entirely offline against recorded archive records, so it
needs no network and no credentials. 122 tests; locally verified on CPython 3.11
and 3.13, both with and without geofetch installed. CI additionally exercises
3.10 and 3.12.

## License

MIT
