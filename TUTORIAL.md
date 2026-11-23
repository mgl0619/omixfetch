# omixfetch tutorial

From a clean machine to downloaded data, raw and processed, with the
metadata alongside.

omixfetch downloads. It doesn't convert between formats. The conversion
commands that still ship are covered once, in the
[appendix](#appendix-transitional-conversion-commands), and are due to move to
a separate tool.

The `omixfetch` outputs shown were captured from live runs; some are trimmed.
The `scanpy` and Seurat snippets in section 3 show how the layout is meant to
be used; they weren't captured, since neither library is a dependency.

---

## 1. Install

**Python 3.10 or newer.** Python 3.9 is past end of life and `pip install`
refuses it. If you see `ModuleNotFoundError: No module named 'omixfetch'`
right after installing, check this first:

```bash
python -c "import sys; print(sys.version)"
```

### With mamba or conda (recommended)

```bash
git clone https://github.com/mgl0619/omixfetch.git
cd omixfetch
mamba env create -f environment.yml      # Python 3.13 + aria2
mamba activate omixfetch
pip install -e ".[dev]"
```

Using `micromamba` for the first time? It needs a root prefix before
`env create` works:

```bash
micromamba shell init -s zsh -p ~/micromamba   # then restart your shell
```

### With a plain virtualenv

```bash
python3.13 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
brew install aria2          # macOS; apt install aria2 on Debian/Ubuntu
```

### Optional transfer engines

```bash
conda install -c conda-forge s5cmd      # S3-native; fast against NCBI's AWS mirror
conda install -c hcc aspera-cli         # UDP transfers from ENA
```

### Check it worked

```bash
omixfetch --version
```

```
omixfetch 0.9.0
aria2c: aria2 version 1.37.0
s5cmd:  not installed
```

If `aria2c` says `not installed`, downloads still work but use a built-in
downloader that can't split a single file across connections. Install aria2
before fetching anything large.

The test suite needs no network:

```bash
python -m pytest
```

---

## 2. Look before you download

Never start a multi-gigabyte transfer blind. `info` summarises a study:

```bash
omixfetch info GSE164378
```

```
 Title         Integrated analysis of multimodal single-cell data
 Organisms     Homo sapiens
 Samples       9
 Files         27
 Resolver      geofetch
 Total size    1.4 GB
 PubMed        34062119
```

`list` shows exactly what would be fetched, and takes the same filters as
`get`:

```bash
omixfetch list GSE164378 --include "*RNA_3P*"
```

`get --dry-run` gives the same preview with the final paths filled in.

---

## 3. Download processed data

For most analyses this is all you need, and it's usually small.

```bash
omixfetch get GSE164378 -o ./data --include "*RNA_3P*" -j 4 -x 16
```

- `-j 4`: four files at once
- `-x 16`: sixteen connections **within** each file

`-x` is the flag that saves the time. Archives throttle each connection, so a
large file arrives much faster split sixteen ways.

Files are grouped by sample, and 10x triplets get their Cell Ranger names:

```
data/GSE164378/
├── data/
│   ├── GSM5008738/
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   └── GSM5008739/
│       └── ...
└── metadata/
    ├── GSE164378.metadata.json
    ├── GSE164378.samples.tsv
    ├── GSE164378.files.tsv
    └── GSE164378.soft.txt
```

`--layout flat` keeps one directory and the published filenames.
`<GSE>_RAW.tar` is skipped when its contents are also listed one by one;
`--include-archives` fetches it anyway.

### Load it

The directory is already in Cell Ranger shape, so there's nothing to rename:

```python
import scanpy as sc
adata = sc.read_10x_mtx("data/GSE164378/data/GSM5008738")
```

```r
library(Seurat)
counts <- Read10X("data/GSE164378/data/GSM5008738")
```

Pair it with the sample annotation:

```python
import pandas as pd
samples = pd.read_csv("data/GSE164378/metadata/GSE164378.samples.tsv", sep="\t")
```

---

## 4. Download raw data

Raw reads come down **as they were deposited**: FASTQ if the submitter
uploaded FASTQ, the original BAM if they uploaded BAM, or NCBI's `.sra`
archive if you ask for it. omixfetch moves and verifies the bytes; it doesn't
change their format.

### 4.1 Preview

Read sets are large. Always preview:

```bash
omixfetch get GSE52778 -o ./data --reads-only --dry-run
```

```
GSE52778  39 file(s), 46.9 GB  ->  ./data/GSE52778
    would fetch  GSM1275863/SRR1039509_1.fastq.gz  (1.1 GB)
    would fetch  GSM1275863/SRR1039509_2.fastq.gz  (1.1 GB)
```

Drop `--dry-run` to fetch. Runs are filed under the GEO sample they belong to.
`--reads` instead of `--reads-only` fetches processed files too.

### 4.2 Choose the form

| You want | Command |
|---|---|
| FASTQ, where the submitter uploaded FASTQ (default) | `--reads-only` |
| the submitter's own files (BAM, CRAM, …) | `--reads-only --prefer submitted` |
| the same, from NCBI's AWS mirror | `--reads-only --source ncbi --sra-flavour original` |
| NCBI's normalised `.sra` | `--reads-only --source ncbi --sra-flavour sra` |

Every file is checked against the archive's published size and md5. A file
that fails is deleted, not left on disk looking finished.

Try one run before committing to a study:

```bash
omixfetch get GSE184708 -o ./data --reads-only --prefer submitted \
    --include "SRR16036138*" --dry-run
```

### 4.3 Check what the submission contains

The form decides what you can do with the data later, so check it first. This
matters most for 10x data submitted as BAM. The cell barcode and UMI live in
BAM tags, and they're lost when NCBI loads the BAM into `.sra`. ENA's FASTQ
is derived the same way and loses them too. For GSE184708, only the original
BAM still has them.

The per-run check, which reads each BAM's header with a 4 MB range request
rather than the whole file, currently runs as part of the planner in the
`fastq` command. Its dry run downloads nothing and shows, for each run, what
was submitted, what the header says, and which files are worth fetching:

```bash
omixfetch fastq GSE184708 -o ./plan --dry-run
```

Before fetching raw data at all, check GEO's processed files. For GSE184708
they're 1.2 GiB against 1.7 TiB of BAMs.

### 4.4 Choose the source and engine

EBI and NCBI both serve the data. Which is faster depends on your network, so
measure from the machine that will do the downloading:

```bash
omixfetch benchmark GSE184708 --sample 100MB
```

```
NCBI S3 (original)        55.0 MiB/s   -> 9.0 h for the full set
NCBI S3 (.sra)            46.6 MiB/s   -> 4.3 h for the full set
ENA HTTPS                  0.9 MiB/s   -> 544.9 h for the full set
```

Those numbers come from one machine on one day, so they **won't** match
yours. Then fetch from the faster source:

```bash
omixfetch get GSE184708 -o ./data --reads-only --source ncbi --sra-flavour original
```

Pick the engine with `--backend`:

| `--backend` | Use for |
|---|---|
| `auto` (default) | `s5cmd` when every file is on S3 and it's installed, otherwise `aria2`, then `python` |
| `aria2` | any HTTPS source; `-j` files × `-x` connections per file |
| `s5cmd` | NCBI's AWS mirror |
| `aspera` | ENA, if your site allows outbound UDP 33001 |
| `python` | built in; whole files only |

```bash
omixfetch get GSE184708 -o ./data --reads-only --source ncbi --sra-flavour original \
    --backend s5cmd -j 8 -x 16
```

### 4.5 Projects with no GEO series

Data deposited straight to SRA, ENA or DDBJ has no GSE. Give omixfetch the
BioProject, study, run or BioSample accession instead:

```bash
omixfetch info PRJNA672658
omixfetch get  PRJNA672658 -o ./data --dry-run
```

```
PRJNA672658  22 file(s), 877.4 GiB  ->  ./data/PRJNA672658
    would fetch  SAMN16565580/SRR12921183_1.fastq.gz  (40.8 GiB)
    would fetch  SAMN16565580/SRR12921183_2.fastq.gz  (33.6 GiB)
```

Runs are grouped under their BioSample accession; `samples.tsv` gives each
one's submitter name (here `Normal_5`, `RSA_5`, …) and characteristics. All
the options in 4.2–4.4 work the same way. To try a single run first:

```bash
omixfetch get SRR12921178 -o ./data
```

### 4.6 When no reads come back

```
GSE164378: ENA serves no public runs for this project. Human data under
controlled access (dbGaP) is the usual reason ...
```

This is an access boundary, not a bug. Restricted human data lives in dbGaP
and needs an approved data access request; no tool can fetch it without one.
The SRA and BioProject accessions are recorded in `metadata.json` for that
request.

---

## 5. Metadata

Every download writes metadata alongside the data. To fetch only metadata:

```bash
omixfetch metadata GSE164378 -o ./data --pep
```

| File | Contents |
|---|---|
| `<ACC>.metadata.json` | the complete parsed record |
| `<ACC>.samples.tsv` | one row per sample, characteristics as columns |
| `<ACC>.files.tsv` | every file, with size and URL |
| `<ACC>.runs.tsv` | sequencing runs mapped to samples, with `--reads` |
| `<ACC>.soft.txt` | the unmodified GEO record |
| `pep/` | PEP project, with `--pep` |

`--pep` writes a `config.yaml` and sample CSV that `peppy`, `looper` and
`eido` read directly.

---

## 6. Resume and manifests

Transfers resume. After an interruption, rerun the same command: partial files
continue from where they stopped, and finished files are skipped after their
size is checked.

To run the transfer with other tooling, export a manifest instead:

```bash
omixfetch manifest GSE164378 -o files.txt --format aria2
aria2c -i files.txt -j 8 -x 16
```

`--format tsv` and `--format json` are also available.

---

## 7. Choosing a resolver

Metadata resolution has two backends, chosen with `--resolver`:

| Value | Behaviour |
|---|---|
| `auto` (default) | geofetch for GSE series, except where it can't serve the request |
| `native` | the built-in SOFT parser, always |
| `geofetch` | geofetch, always (GSE only) |

`auto` falls back to the native parser for GSM and GPL accessions, for any
`--reads` workflow, and for series with no processed files. Those limits were
measured against geofetch 0.12.11, not assumed. geofetch only supplies
metadata here; omixfetch's own engines move every byte.

`--verbose` shows the resolver's own progress logging, which is hidden by
default.

---

## 8. Troubleshooting

**`ModuleNotFoundError: No module named 'omixfetch'`**: the editable install
didn't take effect, usually because `pip` and `pytest` are different
interpreters, or the active Python is 3.9. Use `python -m pytest`, which
always runs the interpreter you think it does.

**`aria2c not found; falling back to the built-in downloader`**: install
aria2. The fallback works, but each file goes over a single connection.

**aria2c is installed but every transfer still falls back**: some networks
put an intercepting TLS proxy in front of the internet, and aria2's TLS
backend refuses it. Check with:

```bash
aria2c --console-log-level=info --max-tries=1 -d /tmp -o t.bin \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164378/suppl/GSE164378_Antibody_HTO_info.xlsx
```

A handshake failure here comes from the environment, not from omixfetch.

**A large download keeps failing**: each retry resumes from the partial file,
so rerunning the same command makes progress.

**ENA is very slow**: run `omixfetch benchmark` and try `--source ncbi`. On
the machines tested so far, NCBI's AWS mirror was roughly 20–60× faster than
ENA's HTTPS.

**Two modalities in one sample kept their long filenames**: this is
intentional. When renaming would make two files collide in one sample
directory, the published names are kept rather than risk overwriting data.

---

## Appendix: transitional conversion commands

These predate the decision to keep omixfetch download-only. They still work,
won't get new features, and will move to a separate conversion tool that reads
omixfetch's manifest.

```bash
omixfetch setup                                    # installs SRA Toolkit + 10x bamtofastq
omixfetch fastq GSE184708 -o ./data                # download, convert, verify, per run
omixfetch get GSE184708 -o ./data --reads-only --source ncbi --convert fastq -t 8
```

`omixfetch fastq` follows the rules in section 4.3. For a 10x BAM submission
it fetches the original BAM and runs `bamtofastq`, and refuses the `.sra`.
Progress is saved after each step, so rerunning skips finished runs. It writes
`download_report.html` and `fastq_manifest.tsv`.

Known issues, which are why conversion is being split out:

- **sra-tools and aria2 can't share a conda environment on Apple Silicon**
  (a `libxml2`/`icu` conflict). `omixfetch setup` installs prebuilt binaries
  outside conda. Alternatively, use `environment-convert.yml` and set
  `OMIXFETCH_FASTERQ_DUMP` to its `fasterq-dump`.
- **`perl-xml-libxml` can't be satisfied**: sra-tools 3.2.x has no osx-arm64
  build. Ask for `sra-tools>=3.4` with `conda-forge` ahead of `bioconda`.
- **`fasterq-dump`: "can't open ... as a RefSeq or as a WGS"**: aligned
  `.sra` files store reads as differences against a reference genome kept in
  NCBI's repository, so conversion needs network access to it. This failed in
  testing on SRR390728. The original BAM is self-contained.
- **`bamtofastq` on Apple Silicon** is an x86_64 binary and needs Rosetta 2
  (`softwareupdate --install-rosetta`).
