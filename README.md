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

For each run it uses the first route that works, falling through on failure:
ENA's own FASTQ, then NCBI's `.sra` converted with `fasterq-dump`, then the
original submission converted with 10x `bamtofastq`.

**10x submissions are handled specially.** When a 10x BAM is loaded into SRA,
the cell barcode and UMI — stored in BAM tags, not as a read — are dropped, so
the `.sra` holds only the cDNA read. GSE184708 shows it plainly: declared
paired-end, 98 bases per spot, `spots_with_mates` 0. Converting that `.sra`
would give FASTQ that Cell Ranger cannot use. The pipeline detects this from
NCBI's metadata, skips the `.sra` route for such runs, and uses the original
BAM with `bamtofastq`, which rebuilds R1/R2/I1 from the tags. Before each BAM
download it reads the header (a few MB) to confirm the `10x_bam_to_fastq`
lines `bamtofastq` needs are present. It checks spot counts
against ENA's metadata, refuses a route whose converter is missing *before*
downloading anything, checks free disk space per run, and records progress in
`<ACC>/.omixfetch/state.json` — rerun the same command after an interruption
and finished runs are skipped. Results are listed in `<ACC>/fastq_manifest.tsv`.

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

| Route | Status |
|---|---|
| 1 GEO processed | real downloads, sizes and gzip integrity checked |
| 2 ENA FASTQ | resolved and range-capable; long transfers were cut off in the test environment |
| 3 ENA submitted | one `.bai` from GSE184708 downloaded, md5 matched ENA |
| 3 via Aspera | endpoints published for all 20 runs of GSE184708; `ascp` not yet run |
| 4 NCBI original | range requests and throughput sampled; no full BAM downloaded |
| 5 NCBI `.sra` | SRR390728 downloaded, md5 matched; its conversion failed on a reference fetch |
| `s5cmd` backend | unit-tested only; not yet run against S3 |
| `fastq` pipeline | unit-tested offline; plans GSE184708 correctly (all 20 runs to `original`); full runs pending |
| bundled tools | SRA Toolkit 3.4.1 and bamtofastq 1.4.1 installed and run on macOS arm64 |
| 10x header probe | GSE184708 BAM header read with a 4 MB range request; `10x_bam_to_fastq` lines found |

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
