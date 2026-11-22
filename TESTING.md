# Testing

## Offline unit tests

```bash
python -m pytest
```

No network, no credentials; every archive response is a recorded fixture.
Live tests are skipped unless asked for.

## Live tests, in cost tiers

Never iterate against a large series. Each tier exists so a change can be
checked at the cheapest level that exercises it.

| Tier | Series | Cost | Exercises |
|---|---|---|---|
| quick | GSE204978, runs SRR19429286 + SRR19429291 | ~1 min, ~10 MB | ENA-FASTQ route; `.sra` route with real `fasterq-dump`; spot counts; cleanup; resume |
| quick | GSE164378 | seconds | controlled-access (dbGaP) series reported, nothing fetched |
| medium | GSE220167, run SRR22543607 | ~30 min, ~25 GB disk | 10x BAM to `bamtofastq`; R1 carries barcode + UMI |
| final | GSE184708 | many hours, 1.7 TiB | the full 20-run 10x series, once, on the target cluster |

```bash
OMIXFETCH_LIVE=quick  python -m pytest tests/live -v
OMIXFETCH_LIVE=medium python -m pytest tests/live -v
OMIXFETCH_LIVE=final  OMIXFETCH_LIVE_OUT=/scratch/me/gse184708 python -m pytest tests/live -v
```

### Why these series

- **GSE204978** has 384 tiny runs with both ENA FASTQ and unaligned `.sra`,
  so both cheap routes can be tested in seconds on the same data.
- **GSE220167** holds the smallest Cell Ranger BAM found among recent 10x
  submissions (8.27 GB), with the `10x_bam_to_fastq` header lines
  `bamtofastq` needs.
- **GSE217662** is a useful counter-example kept in mind rather than tested
  end to end: NCBI labels its originals `TenX`, but the BAM header has no
  `10x_bam_to_fastq` lines, so `bamtofastq` would fail. The header probe
  rejects it before download.
- **GSE164378** is human PBMC data deposited in dbGaP.
- **GSE184708** is the motivating dataset: 20 runs of 10x mouse data
  submitted as BAM, where the `.sra` has lost the barcode read.

## Results

| Tier | Last run | Outcome |
|---|---|---|
| quick | 0.7.x, macOS arm64 | ENA-FASTQ 2/2 in 28 s, sizes match ENA; `.sra` 2/2, spots 23,959 / 26,554 match NCBI in both mates; work dir removed; rerun skips both; GSE164378 reported, exit 1 |
| medium | not yet run | |
| final | not yet run | |
