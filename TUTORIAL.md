# omixfetch tutorial

A walkthrough from a clean machine to a loaded single-cell matrix.

Every `omixfetch` command below has been run against the live archives, and
the outputs shown are what they returned. The `scanpy` and `Seurat` snippets
in section 3 are the intended usage of the directory layout rather than
captured runs — neither library is a dependency of this project.

---

## 1. Install

**Python 3.10 or newer is required.** Python 3.9 is past end of life, and
`pip install` will refuse it — if you see `ModuleNotFoundError: No module
named 'omixfetch'` straight after installing, check this first:

```bash
python -c "import sys; print(sys.version)"
```

### With mamba or conda (recommended)

```bash
git clone https://github.com/mgl0619/omixfetch.git
cd omixfetch
mamba env create -f environment.yml
mamba activate omixfetch
pip install -e ".[dev]"
```

`environment.yml` pins Python 3.13 and includes `aria2`, the download engine.

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

### Check it worked

```bash
omixfetch --version
```

```
omixfetch 0.5.0
aria2c: aria2 version 1.37.0
```

If the second line says `not installed`, downloads still work but fall back to
a slower threaded downloader that cannot split a single file across
connections. Install aria2 before fetching anything large.

Run the test suite to confirm the install is sound — it needs no network:

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

The `Resolver` row tells you which metadata backend ran — see section 7.

`list` shows exactly what would be fetched, and accepts the same filters as
`get`:

```bash
omixfetch list GSE164378 --include "*RNA_3P*"
```

Add `--dry-run` to `get` for the same preview with the final paths resolved.

---

## 3. Download processed data

```bash
omixfetch get GSE164378 -o ./data --include "*RNA_3P*" -j 4 -x 16
```

- `-j 4` — four files at once
- `-x 16` — sixteen connections **within** each file

The second flag is where the time is saved. NCBI throttles per connection, so
a 1 GB matrix downloads far faster split sixteen ways.

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

Pass `--layout flat` to keep one directory with the published filenames.

### Load it

Because the directory is already in Cell Ranger shape, no renaming step is
needed:

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

## 4. Download raw reads

omixfetch fetches FASTQ from **ENA**, which mirrors SRA and serves reads as
ordinary files. No SRA Toolkit, and no `.sra` conversion step.

Always preview first — read sets are large:

```bash
omixfetch get GSE52778 -o ./data --reads-only --dry-run
```

```
GSE52778  39 file(s), 46.9 GB  ->  ./data/GSE52778
    would fetch  GSM1275863/SRR1039509_1.fastq.gz  (1.1 GB)
    would fetch  GSM1275863/SRR1039509_2.fastq.gz  (1.1 GB)
```

Then drop `--dry-run`. Runs are grouped under the GEO sample they belong to.

ENA publishes an md5 for every read file, so these are verified end to end —
a corrupt file is deleted rather than left on disk.

For long-read studies where the submitter deposited BAM rather than FASTQ:

```bash
omixfetch get <GSE> -o ./data --reads-only --prefer submitted
```

### Going faster

Reads are served by both EBI and NCBI, and which is quicker depends on your
network. Measure before committing to a large transfer:

```bash
omixfetch benchmark GSE184708 --sample 100MB
```

```
NCBI S3 (original)        55.0 MiB/s   -> 9.0 h for the full set
NCBI S3 (.sra)            46.6 MiB/s   -> 4.3 h for the full set
ENA HTTPS                  0.9 MiB/s   -> 544.9 h for the full set
```

Those numbers are from one particular machine and are **not** transferable —
run it yourself on the node that will do the downloading.

Then pick a source:

```bash
omixfetch get GSE184708 -o ./data --reads-only --source ncbi
```

If you want FASTQ, the `.sra` route is both the smallest transfer and the
direct path to it — one command downloads and converts:

```bash
omixfetch get GSE184708 -o ./data --reads-only \
    --source ncbi --sra-flavour sra --convert fastq -t 8
```

That needs the SRA Toolkit for the conversion step only; the download itself
does not use it. Install it in its own environment — on Apple Silicon it
cannot coexist with aria2:

```bash
mamba env create -f environment-convert.yml
export OMIXFETCH_FASTERQ_DUMP="$(conda run -n omixfetch-convert which fasterq-dump)"
```

Aspera is worth trying if your site permits outbound UDP 33001:

```bash
conda install -c hcc aspera-cli
omixfetch get GSE184708 -o ./data --reads-only --backend aspera
```

### When no reads come back

```
GSE164378: ENA serves no public runs for this project. Human data under
controlled access (dbGaP) is the usual reason ...
```

This is a permissions boundary, not a bug. Restricted human data lives in
dbGaP and requires an approved data access request; no tool can fetch it
without authorisation. The SRA and BioProject accessions are recorded in
`metadata.json` so you can take them to `prefetch` once you have access.

---

## 5. Metadata

Every download writes metadata alongside the data. To fetch only that:

```bash
omixfetch metadata GSE164378 -o ./data --pep
```

| File | Contents |
|---|---|
| `<ACC>.metadata.json` | the complete parsed record |
| `<ACC>.samples.tsv` | one row per sample, characteristics as columns |
| `<ACC>.files.tsv` | file manifest with sizes and URLs |
| `<ACC>.runs.tsv` | sequencing runs, with `--reads` |
| `<ACC>.soft.txt` | the unmodified GEO record |
| `pep/` | PEP project, with `--pep` |

`--pep` writes a `config.yaml` plus sample CSV that `peppy`, `looper` and
`eido` read directly.

---

## 6. Resume and manifests

Transfers resume. Rerun the same command after an interruption and aria2
continues from where it stopped — partial files are never silently treated as
complete.

To hand the transfer to another system, emit a manifest instead:

```bash
omixfetch manifest GSE164378 -o files.txt --format aria2
aria2c -i files.txt -j 8 -x 16
```

`--format tsv` and `--format json` are also available.

---

## 7. Choosing a resolver

Metadata resolution has two backends, selected with `--resolver`:

| Value | Behaviour |
|---|---|
| `auto` (default) | geofetch for GSE series, except where it cannot serve the request |
| `native` | the built-in SOFT parser, always |
| `geofetch` | geofetch, always (GSE only) |

`auto` falls back to native for GSM and GPL accessions, for any `--reads`
workflow, and for series with no processed files. Those limits were measured
against geofetch 0.12.11, not assumed. Transfers go through aria2c either way:
geofetch is called in metadata-only mode and never downloads anything.

Use `--verbose` to see the resolver's own progress logging, which is
suppressed by default.

---

## 8. Troubleshooting

**`ModuleNotFoundError: No module named 'omixfetch'`** — the editable install
did not take effect, usually because `pip` and `pytest` are different
interpreters, or the active Python is 3.9. Use `python -m pytest`, which
always uses the interpreter you think it does.

**`aria2c not found; falling back to the built-in downloader`** — install
aria2. The fallback works but transfers each file over a single connection.

**aria2c installed but every transfer still falls back** — some environments
put an intercepting TLS proxy in front of the network that aria2's TLS backend
refuses. Verify with:

```bash
aria2c --console-log-level=info --max-tries=1 -d /tmp -o t.bin \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164378/suppl/GSE164378_Antibody_HTO_info.xlsx
```

A handshake failure here is environmental, not an omixfetch bug.

**A large download fails repeatedly** — each retry resumes from the partial
file, so simply rerunning the same command makes forward progress.

**`sra-tools` will not install / `perl-xml-libxml` cannot be satisfied** —
on osx-arm64, sra-tools 3.2.x depends on a package with no Apple Silicon
build. Ask for a newer one (`sra-tools>=3.4`) and put `conda-forge` ahead of
`bioconda`.

Separately, sra-tools and aria2 cannot be solved into one environment on that
platform at all — a combined solve fails on a `libxml2`/`icu` conflict. Use
`environment-convert.yml` to keep them apart. If you have previously installed
sra-tools into an environment that had aria2, check `aria2c --version` still
works there: an environment seen here ended up with sra-tools present and
aria2 gone, though which install removed it was not established.

**`fasterq-dump` fails with "can't open ... as a RefSeq or as a WGS"** — the
`.sra` does not contain the read bases. Aligned submissions are stored as
differences against a reference genome, the same idea as CRAM, and the
reference lives in NCBI's repository rather than in the file you downloaded.
Installing the toolkit locally does not help: the software is present, the
data is incomplete. This is also where part of the size advantage over the
original BAM comes from. Verified failing this way on SRR390728.

Three ways out, in order of preference:

1. **Let it fetch, and cache sensibly.** References are cached in the VDB user
   repository, so one reference set serves every run of a study — a few GB
   once, not per run. Point the cache at scratch, never at a home directory
   with a quota, or a failed cache write will look like a network failure:

   ```bash
   vdb-config --interactive     # set the cache location
   ```

2. **Pre-fetch.** `prefetch SRR16036138` acquires the references alongside the
   run into the local repository, after which conversion works offline.

3. **Avoid references entirely** with `--sra-flavour original`: submitted
   BAMs are self-contained. For GSE184708 that is 1.7 TiB instead of 709 GiB,
   and 10x BAMs then need `bamtofastq` rather than `fasterq-dump`.

The caching behaviour in (1) is how the toolkit is documented to work; it was
not verified here, because the test never got past the first reference fetch.

**Two modalities in one sample kept their long filenames** — intentional. When
renaming would make two files collide inside one sample directory, the
published names are kept rather than risk overwriting data.
