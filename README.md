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

## Why this exists

Downloading a GEO series is nominally easy and repeatedly annoying:

- GEO's HTTPS endpoint returns **403 for directory listings**, so you cannot
  simply enumerate `.../suppl/`. The authoritative file list lives in the SOFT
  record.
- NCBI throttles bursts by returning an **HTML error page with a 200 status**,
  which naive clients happily write to disk as a `.tar`.
- A 40 GB `_RAW.tar` over one connection is slow; GEO supports range requests,
  so it does not have to be.
- Metadata and data are fetched by different tools, so the characteristics
  table that makes the matrix interpretable often never gets saved at all.
- Getting raw reads usually means a second toolchain (SRA Toolkit) and an
  `.sra` → FASTQ conversion step, even though ENA mirrors the same reads as
  plain, range-capable HTTPS files.
- `<GSE>_RAW.tar` is a bundle of the same per-sample supplementary files GEO
  also serves individually. Fetch both and you transfer the data twice — for
  GSE164378 that is 2.94 GB for 1.47 GB of unique data.

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

By default omixfetch parses GEO's SOFT records itself. With the optional extra
it can instead delegate metadata resolution to
[geofetch](https://pep.databio.org/geofetch/):

```bash
pip install "omixfetch[geofetch]"
omixfetch get GSE164378 -o ./data --resolver geofetch
```

Both resolvers agree on GSE164378 (27 files, 1.47 GB). geofetch brings its own
maintained parsing and the pepkit ecosystem; the native resolver keeps the
install at four dependencies. Transfers go through aria2c either way.

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

## What you get

```
data/GSE164378/
├── data/
│   ├── GSE164378_RAW.tar
│   └── GSE164378_sc.meta.data_3P.csv.gz
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

## Development

```bash
mamba env create -f environment.yml && mamba activate omixfetch
pip install -e ".[dev]"
python -m pytest
ruff check src tests
```

The test suite runs entirely offline against recorded archive records, so it
needs no network and no credentials. 82 tests; locally verified on CPython 3.11
and 3.13, both with and without geofetch installed. CI additionally exercises
3.10 and 3.12.

## License

MIT
