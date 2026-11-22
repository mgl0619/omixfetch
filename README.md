# omixfetch

**Fast, verified downloads of public genomics data, with the metadata that
makes it usable.**

Give omixfetch a GEO accession and it works out every file behind it:
processed matrices, raw reads in whatever form they were deposited, and the
study and sample annotation. It picks the fastest source from your machine,
moves the bytes in parallel, checks them against the archive's own sizes and
checksums, and lays them out ready for analysis.

```bash
omixfetch info GSE184708                 # what is this study, and how big?
omixfetch list GSE184708                 # which files would I get?
omixfetch get  GSE184708 -o ./data       # fetch processed files + metadata
```

omixfetch **downloads**. It doesn't convert between formats; see
[Scope](#scope) for how conversion is handled.

---

## Contents

- [What it does](#what-it-does)
- [Install](#install)
- [Downloading processed data](#downloading-processed-data)
- [Downloading raw data](#downloading-raw-data)
- [Know what you are downloading](#know-what-you-are-downloading)
- [Speed: sources and transfer engines](#speed-sources-and-transfer-engines)
- [What you get on disk](#what-you-get-on-disk)
- [Metadata](#metadata)
- [What has been verified](#what-has-been-verified)
- [Scope](#scope)
- [Python API](#python-api) · [Related tools](#related-tools) · [Development](#development)

---

## What it does

| | |
|---|---|
| **Resolves** | a GEO series, sample or platform into its complete file list: processed supplementary files, series matrices, and raw reads from ENA and NCBI |
| **Chooses sources** | the same raw data is served by EBI (ENA) and by NCBI's AWS mirror; omixfetch measures both from your machine and uses the faster |
| **Transfers** | in parallel: many files at once, and many connections per file (`aria2c`), S3-native (`s5cmd`) or UDP (`ascp`) |
| **Verifies** | every file against the archive's published size, and md5 where one is published (ENA and NCBI publish md5 for raw files) |
| **Resumes** | interrupted transfers continue where they stopped; rerunning a finished job fetches nothing |
| **Lays out** | files by sample, with 10x matrices given the names `scanpy` and Seurat expect |
| **Records** | study and sample metadata as JSON, TSV and PEP, plus the unmodified archive records |
| **Avoids waste** | drops `_RAW.tar` when its contents are listed individually, checks a BAM's header before committing to a 100 GB download, and reports controlled-access data instead of failing on it |

Supported today: **NCBI GEO**, with raw reads via **ENA** and **NCBI SRA**.
Planned: ArrayExpress/E-MTAB, PRIDE, MetaboLights (see [ROADMAP.md](ROADMAP.md)).

---

## Install

Python 3.10 or newer. Python 3.9 is past end of life and `pip` will refuse it.

```bash
git clone https://github.com/mgl0619/omixfetch.git && cd omixfetch
mamba env create -f environment.yml      # or micromamba / conda; includes aria2
mamba activate omixfetch
pip install -e .
omixfetch --version
```

`aria2c` does the parallel transfers and is strongly recommended. Without it,
omixfetch falls back to a built-in downloader that fetches several files at
once but can't split a single file across connections. Plain virtualenv
instead of conda:

```bash
python3.13 -m venv .venv && source .venv/bin/activate
pip install -e .
brew install aria2                       # macOS;  apt install aria2  on Debian/Ubuntu
```

Optional engines: `conda install -c conda-forge s5cmd` (fast for NCBI's AWS
mirror), `conda install -c hcc aspera-cli` (ENA over UDP).

> `ModuleNotFoundError: No module named 'omixfetch'` straight after installing
> almost always means `pip` and `pytest` are different interpreters, or the
> active Python is 3.9. Use `python -m pytest`, which always runs the
> interpreter you think it does.

---

## Downloading processed data

This is what most analyses need, and it's usually small. GSE184708, a 20-run
10x study with 1.7 TiB of raw BAMs, ships its count matrix on GEO at 1.2 GiB.

```bash
omixfetch info GSE184708                          # title, organism, samples, total size
omixfetch list GSE184708                          # every file, with sizes
omixfetch get  GSE184708 -o ./data                # fetch them, plus metadata
```

Filter what you fetch:

```bash
omixfetch get GSE164378 -o ./data \
    --include "*.mtx.gz" --include "*.tsv.gz" \  # globs; repeatable
    --exclude "*_RAW.tar" \
    --max-size 5GB \
    --matrix                                      # also the series-matrix file
omixfetch get GSE164378 -o ./data --dry-run       # show what would be fetched
omixfetch metadata GSE164378 -o ./data            # metadata only
```

`<GSE>_RAW.tar` is dropped automatically when its member files are also listed
individually: it's a bundle of the same data, and fetching both would move
everything twice. `--include-archives` keeps it.

---

## Downloading raw data

Raw reads are fetched **as they were deposited**: FASTQ if the submitter
uploaded FASTQ, the original BAM if they uploaded BAM, or NCBI's `.sra`
archive. Add `--reads` to fetch them alongside processed files, or
`--reads-only` to fetch only them.

```bash
omixfetch get GSE52778 -o ./data --reads-only --dry-run    # always preview first
omixfetch get GSE52778 -o ./data --reads-only
```

### Which raw files exist

Each run is served in up to three forms. Here is SRR16036138 from GSE184708:

| Form | What it is | Size | Source | Select with |
|---|---|---|---|---|
| **FASTQ** | reads as files, when the submitter uploaded FASTQ | — (this run was submitted as BAM) | ENA | `--source ena` (default) |
| **Original submission** | the exact file uploaded, here a 10x BAM | 119.2 GB | ENA, NCBI AWS | `--prefer submitted`, or `--source ncbi --sra-flavour original` |
| **Normalised `.sra`** | NCBI's archive format | 48.1 GB | NCBI AWS | `--source ncbi --sra-flavour sra` |

NCBI also offers **SRA Lite** (22.2 GB here), the same reads with quality
scores reduced to pass/fail. omixfetch doesn't fetch it yet; see the roadmap.

```bash
# FASTQ from ENA, when it exists
omixfetch get GSE52778  -o ./data --reads-only

# the original submission, from ENA or from NCBI's AWS mirror
omixfetch get GSE184708 -o ./data --reads-only --prefer submitted
omixfetch get GSE184708 -o ./data --reads-only --source ncbi --sra-flavour original

# NCBI's .sra
omixfetch get GSE184708 -o ./data --reads-only --source ncbi --sra-flavour sra
```

Add `--include 'SRR16036138*'` to try a single run before committing to a
whole study.

### Controlled access

Human data deposited in dbGaP has a BioProject but no public reads.
omixfetch says so rather than reporting an empty download:

```
GSE164378: ENA serves no public runs for this project. Human data under
controlled access (dbGaP) is the usual reason ...
```

No tool can fetch those without an approved data access request. The SRA and
BioProject accessions are recorded in the metadata for that request.

---

## Know what you are downloading

Raw downloads are large, and the form you fetch decides what you can do with
it later. Two things omixfetch checks before you spend hours on a transfer:

**1. Was the submission FASTQ or BAM?** omixfetch reads NCBI's record of the
original file, backed by ENA's. For a FASTQ submission, every form contains
every read. For a BAM submission, the `.sra` and ENA's FASTQ are *derived from
the BAM* and carry whatever it lost.

**2. For a BAM, what is in it?** omixfetch reads the BAM header with a 4 MB
range request, not the whole file, and reports whether the BAM still holds
every sequenced read:

| Header shows | Means |
|---|---|
| `10x_bam_to_fastq` lines | a Cell Ranger BAM; barcodes and UMIs are recoverable from its tags |
| aligner kept unmapped reads (STAR `--outSAMunmapped Within`, BWA, minimap2, Bowtie2/HISAT2 without `--no-unal`) | complete |
| aligner dropped unmapped reads | contains only what mapped |
| single-cell, but no 10x lines | barcodes cannot be recovered as reads |

**The 10x case, which matters most.** When a 10x BAM is loaded into SRA, the
cell barcode and UMI, stored as BAM tags rather than as a read, are dropped.
GSE184708 shows it plainly in NCBI's RunInfo for SRR16036138:
`LibraryLayout` PAIRED, `spots_with_mates` **0**, `avgLength` **98**, which is
the cDNA read alone. Its `.sra` and ENA's FASTQ therefore **cannot be used for
single-cell analysis**. Only the original BAM still has the barcodes. Its
header confirms they are there (`10x_bam_to_fastq:R1(CR:CY,UR:UY)`) and that
unmapped reads were kept (`--outSAMunmapped Within`).

NCBI's own labels aren't enough to tell. GSE184708 is labelled
`TRANSCRIPTOMIC`, not single-cell, even though it is 10x data. GSE217662 is
labelled `TenX`, but its BAM headers have no rebuild lines. That's why the
header itself is read.

**Before fetching raw data at all, check GEO's processed files.** For
GSE184708 they are 1.2 GiB against 1.7 TiB. Raw data is worth fetching for
reprocessing (a new reference, a different aligner, ambient-RNA correction)
or for tools that read alignments directly (velocyto, genotype
demultiplexing), not for re-analysing the authors' counts.

---

## Speed: sources and transfer engines

Which source is faster depends entirely on your network. Measure it from the
machine that will do the downloading:

```bash
omixfetch benchmark GSE184708 --sample 100MB
```

It samples a fixed byte range from each source and reports MiB/s plus the
projected time for the whole study, costing a few hundred MB rather than
terabytes. One measurement from development: ENA 0.5 MiB/s against NCBI's AWS
mirror 11.9 MiB/s. Yours will differ, so measure.

| `--backend` | Use for | Parallelism |
|---|---|---|
| `aria2` | any HTTPS source | files × connections per file (`-j`, `-x`) |
| `s5cmd` | NCBI's AWS mirror | objects × multipart parts |
| `aspera` | ENA; needs outbound UDP 33001 | UDP transfer |
| `python` | anywhere; built in | files only |

`auto` (the default) picks `s5cmd` when every file is on S3 and it is
installed, otherwise `aria2`, falling back to `python` if `aria2c` fails.

```bash
omixfetch get GSE184708 -o ./data --reads-only --source ncbi --sra-flavour original \
    --backend s5cmd -j 8 -x 16
```

To run the transfer with your own tooling, export a manifest:

```bash
omixfetch manifest GSE164378 -o files.txt --format aria2   # or tsv, json
aria2c -i files.txt -j 8 -x 16
```

---

## What you get on disk

Files are grouped by sample. 10x matrix triplets get their Cell Ranger names,
so `scanpy.read_10x_mtx()` and Seurat's `Read10X()` read a directory directly:

```
data/GSE164378/
├── data/
│   ├── GSM5008738/
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── GSM5008739/
│   │   └── ...
│   └── GSE164378_sc.meta.data_3P.csv.gz     series-level files stay at the top
└── metadata/
    ├── GSE164378.metadata.json              complete parsed record
    ├── GSE164378.samples.tsv                one row per sample
    ├── GSE164378.files.tsv                  every file, with size and URL
    ├── GSE164378.runs.tsv                   sequencing runs (with --reads)
    ├── GSE164378.soft.txt                   unmodified GEO records
    ├── GSE164378.samples.soft.txt
    └── pep/                                 with --pep
        ├── GSE164378_config.yaml
        └── GSE164378_samples.csv
```

`--layout flat` keeps one directory and the published names. Renaming happens
only when unambiguous: if two files in one sample would collide, the original
names are kept.

---

## Metadata

`samples.tsv` flattens each sample's characteristics into columns, taking the
union across samples so sparsely annotated samples get empty cells rather
than being dropped. `runs.tsv` maps each sequencing run to its GEO sample,
with read counts and md5s. The unmodified SOFT records are kept, so nothing is
lost to the parser.

`--pep` writes a [PEP](https://pep.databio.org/) project (a `config.yaml` and a
sample CSV) that `peppy`, `looper` and `eido` read directly.

**Resolvers.** GEO metadata is parsed either by omixfetch's own SOFT parser or
by [geofetch](https://pep.databio.org/geofetch/). `--resolver auto` (the
default) uses geofetch for GSE series and the built-in parser for GSM and GPL
accessions, which geofetch's project API does not handle, and for read
workflows. `omixfetch info` shows which one ran.

---

## What has been verified

Run against the live archives, not mocked. 203 offline unit tests cover the
logic; [TESTING.md](TESTING.md) describes the live test tiers.

| | Status | Evidence |
|---|---|---|
| GEO processed files | ✓ | real downloads; sizes and gzip integrity checked |
| ENA FASTQ | ✓ | GSE204978, 2 runs; sizes match ENA to the byte, md5 checked |
| NCBI `.sra` from AWS | ✓ | SRR390728 and GSE204978 runs; md5 matched |
| ENA submitted file | ✓ | a `.bai` from GSE184708; md5 matched |
| Source measurement | ✓ | GSE204978: ENA 0.5 vs NCBI 11.9 MiB/s; NCBI chosen |
| BAM header check | ✓ | GSE184708 header read with a 4 MB range request |
| Resume | ✓ | rerun fetched nothing |
| Controlled access | ✓ | GSE164378 reported, exit code 1 |
| Full original BAM from AWS | not yet | range requests and throughput sampled only |
| `aria2c` | on users' machines | the development sandbox blocks its TLS; works on real machines |
| `s5cmd`, Aspera | not yet | unit-tested / endpoints confirmed only |

---

## Scope

omixfetch downloads. Converting between formats, such as `.sra` or BAM to
FASTQ, or `.sra` to BAM, belongs in a **separate tool**. Conversion needs its
own heavy dependencies (SRA Toolkit, 10x `bamtofastq`, `samtools`), has
platform quirks (sra-tools and aria2 can't share a conda environment on Apple
Silicon; aligned `.sra` needs NCBI's reference service), and fails in
different ways from downloading. The planned conversion tool will read
omixfetch's file manifest as its input. See [ROADMAP.md](ROADMAP.md).

**Transitional.** Conversion commands written before this decision still
ship, and keep working until the split:

```bash
omixfetch setup                                   # installs SRA Toolkit + bamtofastq
omixfetch fastq GSE184708 -o ./data --dry-run     # download-and-convert pipeline, as a plan
omixfetch get GSE184708 -o ./data --reads-only --source ncbi --convert fastq
```

They follow the same rules as above: for a 10x BAM submission they use the
original BAM with `bamtofastq` and refuse the `.sra`, since it has lost the
barcodes. They write an HTML report (`download_report.html`) of what was
decided and why; examples are in [docs/examples](docs/examples). These
commands won't receive new features.

---

## Python API

```python
from omixfetch import resolve
from omixfetch.backends import Aria2Backend

dataset = resolve("GSE164378")
print(dataset.metadata.title, dataset.metadata.n_samples, dataset.total_bytes)

big = [f for f in dataset.files if (f.size or 0) > 1e9]
Aria2Backend(connections=16, concurrent=4).download(big, "./data")
```

---

## Related tools

[geofetch](https://pep.databio.org/geofetch/) (GEO + SRA, PEP output; used by
omixfetch as a resolver), [GEOquery](https://bioconductor.org/packages/GEOquery/)
(R), [GEOparse](https://geoparse.readthedocs.io/) (Python SOFT parsing),
[ffq](https://github.com/pachterlab/ffq) (links and metadata only) and
[fastq-dl](https://github.com/rpetit3/fastq-dl) (ENA/SRA reads). omixfetch's
focus is the download itself: choosing the source, moving the bytes in
parallel, verifying them, and knowing before a 100 GB transfer whether the
file is the one you need.

---

## Development

```bash
pip install -e ".[dev]"
python -m pytest                                   # offline, ~1 s
OMIXFETCH_LIVE=quick python -m pytest tests/live    # real archives, ~2 min
ruff check src tests
```

[TUTORIAL.md](TUTORIAL.md) walks through a full session.
[TESTING.md](TESTING.md) lists the live test tiers and results.
[ROADMAP.md](ROADMAP.md) has planned archives and deferred work.

## License

MIT
