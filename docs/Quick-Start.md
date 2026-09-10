# Quick start

## 0. Never run on a login node

Every command below goes through SLURM. `bin/run_pipeline.sbatch` runs the
Nextflow *driver* inside a small allocation; the driver then submits each task
as its own job. Everything after the script name is passed to `nextflow run`.

```bash
sbatch bin/run_pipeline.sbatch -profile uwa,apptainer --fasta ... --outdir ...
```

On a busy cluster, override the wall time so the driver backfills rather than
waiting for a priority slot — see
[Troubleshooting](Troubleshooting.md#jobs-stuck-pending-priority).

## 1. Check your setup

```bash
nextflow -version          # need >= 24.04
apptainer --version        # or singularity / docker / conda
sinfo -p work -o "%.6D %.20C"
```

## 2. Run the built-in test

The `test` profile runs every stage on a reference and read set small enough
to finish in minutes. It is the fastest way to confirm your container runtime,
filesystem bindings and scheduler settings work before committing to a real run.

```bash
python3 tests/make_test_data.py -o tests/data     # generates a synthetic cohort

sbatch --time=04:00:00 --cpus-per-task=1 --mem=4G \
    bin/run_pipeline.sbatch -profile test,apptainer \
    --input tests/test_samplesheet.csv \
    --fasta tests/data/test_ref.fna \
    --outdir test_results
```

The generated cohort gives `SAMPLE_A` two sequencing runs, so the test covers
the per-sample merge path as well as everything else.

You should see `results/variants/test.filtered.vcf.gz` and a populated
`test_results/benchmarks/` directory.

## 3. Fetch the reads

The full cohort is ~0.71 TB across 476 files. Two options.

**As a chain of short jobs** — recommended on a busy cluster. Each link
resumes where the last stopped, and the chain exits early once every file
verifies:

```bash
bin/fetch_reads_chain.sh /group/peg/cicer/cret/reads 12 3 8
#                        <outdir>                     n  h  parallel
```

Twelve 1-core, 3-hour jobs. Small and short enough to backfill, where a single
long job waits for a priority slot.

**Inside the pipeline** — anything missing from `--reads_dir` is downloaded and
checksum-verified, and `storeDir` means a re-run never re-downloads. Convenient,
but it ties the download to the pipeline's own wall time.

**Verify or stop at any point:**

```bash
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
for n in 2 8 24; do
  sbatch bin/run_pipeline.sbatch -profile uwa,apptainer,bench_${n} \
      --sra_metadata assets/cret_metadata.tsv \
      --fasta /group/peg/cicer/chickpea/genome/PBA_HatTrick/PBA_HatTrick.fasta \
      --reads_dir /group/peg/cicer/cret/reads \
      --outdir results_bench${n} --benchmark_label bench${n}
done
```

## 5. The full run

```bash
sbatch bin/run_pipeline.sbatch -profile uwa,apptainer \
    --sra_metadata assets/cret_metadata.tsv \
    --fasta /group/peg/cicer/chickpea/genome/PBA_HatTrick/PBA_HatTrick.fasta \
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
