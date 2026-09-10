# Quick start

## 1. Check your setup

```bash
nextflow -version          # need >= 24.04
apptainer --version        # or singularity / docker / conda
```

## 2. Run the built-in test

The `test` profile runs every stage on a reference and read set small enough
to finish in minutes. It is the fastest way to confirm your container runtime,
filesystem bindings and scheduler settings work before committing to a real run.

```bash
nextflow run . -profile test,apptainer \
    --input tests/test_samplesheet.csv \
    --fasta tests/data/test_ref.fna \
    --outdir test_results
```

You should see `results/variants/test.filtered.vcf.gz` and a populated
`test_results/benchmarks/` directory.

## 3. Fetch the reads

The full cohort is ~0.71 TB across 476 files. Two options.

**Inside the pipeline** — the default. Anything missing from `--reads_dir`
is downloaded and checksum-verified, and `storeDir` means a re-run never
re-downloads:

```bash
nextflow run . -profile uwa,apptainer \
    --sra_metadata assets/cret_metadata.tsv \
    --fasta /path/to/assembly.fna \
    --reads_dir /group/peg/cicer/cret/reads \
    --outdir results
```

**Outside the pipeline** — useful when you want the download running while you
work on something else. Resumable and idempotent:

```bash
bin/fetch_reads.sh -o /group/peg/cicer/cret/reads -j 8
bin/fetch_reads.sh -o /group/peg/cicer/cret/reads -c   # verify only
bin/stop_fetch.sh  /group/peg/cicer/cret/reads         # stop cleanly
```

> Do not run two fetchers against the same directory. The script takes an
> `flock` and will refuse, but killing one incorrectly can leave orphaned
> `curl` children — always stop with `bin/stop_fetch.sh`, never a bare `kill`.

ENA throttles above roughly 8 concurrent connections; `-j 16` measured
*slower* than `-j 8` in practice. Expect 12–15 hours for the full set.

## 4. Scale up gradually

Do not go straight to 238 runs. Measure first — see
[Benchmarking](Benchmarking.md).

```bash
# 2 samples, then 8, then 24
nextflow run . -profile uwa,apptainer,bench_2 \
    --sra_metadata assets/cret_metadata.tsv \
    --fasta /path/to/assembly.fna \
    --reads_dir /group/peg/cicer/cret/reads \
    --outdir results_bench2 --benchmark_label bench2
```

## 5. The full run

```bash
nextflow run . -profile uwa,apptainer \
    --sra_metadata assets/cret_metadata.tsv \
    --fasta /path/to/assembly.fna \
    --reads_dir /group/peg/cicer/cret/reads \
    --outdir results \
    --benchmark_label full_v1 \
    --run_winpca \
    -resume
```

`-resume` is safe and should be your default. Nextflow caches on input
content, so re-running after a failure recomputes only what actually changed.

## Using a samplesheet instead

If your reads did not come from SRA, skip `--sra_metadata` and supply
`--input` with columns `sample,run,fastq_1,fastq_2`:

```csv
sample,run,fastq_1,fastq_2
CudiB_009,SRR6242246,/path/SRR6242246_1.fastq.gz,/path/SRR6242246_2.fastq.gz
CudiB_009,SRR6242299,/path/SRR6242299_1.fastq.gz,/path/SRR6242299_2.fastq.gz
```

Two rows sharing a `sample` are two runs of that sample and will be merged.
