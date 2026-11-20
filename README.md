# omixfetch

Parallel, resumable downloads of public genomics datasets from **NCBI GEO**, with
full study and sample metadata extracted alongside the data.

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

## Install

```bash
pip install omixfetch
```

`aria2c` is strongly recommended and used automatically when present:

```bash
brew install aria2                        # macOS
sudo apt install aria2                    # Debian/Ubuntu
conda install -c conda-forge aria2        # conda
```

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

**0.1 — GEO.** Series (GSE), samples (GSM) and platforms (GPL): supplementary
files, series matrix files, study metadata, per-sample characteristics.

**0.2 — ArrayExpress.** E-* accessions via the EBI BioStudies API, with
per-sample annotation parsed from MAGE-TAB SDRF files. Work in progress on the
`arrayexpress-support` branch.

Not in scope: downloading raw reads from SRA/ENA. Use
[SRA Toolkit](https://github.com/ncbi/sra-tools) or
[fastq-dl](https://github.com/rpetit3/fastq-dl); omixfetch records the SRA and
BioProject identifiers in its metadata so you can hand them straight over.

## Related tools

[GEOfetch](https://pep.databio.org/geofetch/) (GEO + SRA, PEP metadata output),
[GEOquery](https://bioconductor.org/packages/GEOquery/) (R, the most complete
GEO client), [GEOparse](https://geoparse.readthedocs.io/) (Python, SOFT
parsing) and [ffq](https://github.com/pachterlab/ffq) (metadata and links only).
omixfetch differs in treating the transfer itself as the problem worth solving.

## Development

```bash
pip install -e ".[dev]"
pytest
ruff check src tests
```

The test suite runs entirely offline against recorded archive records.

## License

MIT
