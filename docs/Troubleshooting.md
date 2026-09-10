# Troubleshooting

## Never run on the login node

Every compute step must go through SLURM.

```bash
sbatch bin/run_pipeline.sbatch -profile uwa,apptainer ...
sbatch bin/fetch_reads.sbatch /group/peg/cicer/cret/reads 8
```

`bin/run_pipeline.sbatch` runs the Nextflow *driver* in a small allocation;
the driver submits each task as its own job. The `test`, `uwa` and `setonix`
profiles all set `executor = 'slurm'`, so a plain `nextflow run` with one of
those will still submit tasks correctly — but the driver itself would be on
the login node, so use the sbatch wrapper.

To check nothing of yours is running locally:

```bash
ps -u $USER -o pcpu=,args= --sort=-pcpu | head
```

## Downloads

**Two fetchers at once corrupt files.** `fetch_reads.sh` takes an `flock` and
refuses to start alongside another, but a badly-killed run leaves orphaned
`curl` children that keep writing. Always stop with:

```bash
bin/stop_fetch.sh /group/peg/cicer/cret/reads
```

Never `kill` the script directly — it runs under `setsid` as its own process
group precisely so `stop_fetch.sh` can signal the whole tree. Killing only the
parent leaves the `xargs` layer alive and respawning.

**Verify what you have:**

```bash
bin/fetch_reads.sh -o /group/peg/cicer/cret/reads -c
```

Reports `OK-CACHED`, `BAD-MD5` or `MISSING` per file. Anything failing MD5 is
deleted and refetched automatically on the next real run — truncated files are
expected after an interruption and are not a cause for concern.

**Downloads are slow.** Expect 12–15 hours for 0.71 TB. More parallelism does
*not* help: measured throughput at `-j 16` was about half that at `-j 8`,
because ENA throttles per-connection above roughly 8. Stay at 6–8.

## `A process input channel evaluates to null`

A `val` input received `null`. Nextflow rejects this. Optional parameters must
be given a concrete default before being passed to a process — for example
`chr_regex ?: '.'`. If you add a parameter that is `null` by default and pass
it to a process as `val`, normalise it in the subworkflow first.

## `Failed to pull singularity image` / `manifest unknown`

The image name resolved to Docker Hub but lives on quay.io. All bioconda
containers need the explicit registry:

```groovy
container "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
```

Bare `biocontainers/...` resolves to `docker.io/biocontainers/...`, which for
most of these tags does not exist. `broadinstitute/gatk` and `python` are
genuine Docker Hub images and correctly have no prefix.

If a pull times out, the GATK image is large — raise
`apptainer.pullTimeout`, or pre-pull once into the shared cache:

```bash
export NXF_APPTAINER_CACHEDIR=/group/sae001/ahockey/apptainer_cache
apptainer pull docker://broadinstitute/gatk:4.6.1.0
```

## Out-of-memory kills

`conf/base.config` retries exit codes 104, 134, 137, 139, 140, 143 and 247 up
to three times, scaling memory with `task.attempt`. A task that fails all
three attempts is genuinely under-provisioned — do not simply raise
`maxRetries`. Find the process in the trace:

```bash
bin/summarise_benchmark.py results/benchmarks/trace-*.txt --include-failed
```

and raise that process's `withName` block. `GENOMICSDBIMPORT` and
`GENOTYPEGVCFS` are the usual culprits: their memory grows with cohort size,
so a limit that held at 24 samples may not hold at 161.

## Too many tasks

An assembly with many unplaced scaffolds scatters into one task per sequence.
The *C. echinospermum* assembly has 17,304 sequences, which at 161 samples
would be ~2.8 million `HaplotypeCaller` tasks. Restrict the scatter:

```bash
--chr_regex '^cicec\.S2Drd065\.gnm1\.chr' --intervals_min_length 100000
```

`BUILD_INTERVALS` fails loudly if the regex matches nothing, rather than
silently producing an empty cohort.

## Jobs stuck `PENDING (Priority)`

The cluster is busy, not broken. Check with:

```bash
sinfo -p work -o "%.10P %.6D %.20C"     # CPUS(A/I/O/T): allocated/idle/other/total
squeue -h -p work -t PENDING | wc -l    # how many jobs are ahead of you
```

**`squeue --start` estimates ignore backfill.** A job can show a start time two
days out and then begin within a minute, because the scheduler fits small,
short jobs into gaps ahead of large reservations. Do not resubmit on the basis
of that estimate alone.

**What actually determines whether you backfill** is how big and how long the
job is. A 4-core, 2-day request has to wait for a full priority slot. A
1-core, 3-hour request slots into almost any gap.

This matters most for the download, which is network-bound and needs almost no
CPU. Rather than one long job, submit a chain of short ones:

```bash
bin/fetch_reads_chain.sh /group/peg/cicer/cret/reads 12 3 8
#                        <outdir>                    n  h  parallel
```

Twelve 1-core, 3-hour jobs chained with `--dependency=afterany`. Each resumes
exactly where the last stopped, because `fetch_reads.sh` is resumable and
checksum-verified, and once every file verifies the remaining links exit
within seconds. On a cluster with 187 jobs queued ahead, the first link
started in under a minute where a 4-core/2-day job was estimated 32 hours out.

The same principle applies to the Nextflow driver. Override the wall time at
submission when you expect a short run:

```bash
sbatch --time=04:00:00 --cpus-per-task=1 --mem=4G bin/run_pipeline.sbatch ...
```

Keep `queueSize` and `submitRateLimit` in the profile as they are — they stop
Nextflow from flooding the scheduler with thousands of task submissions.

## `-resume` did not reuse anything

Nextflow caches on input *content*, not filenames. Common causes:

- The work directory changed. Keep `-w` stable across runs.
- A file was re-downloaded and now has a different timestamp *and* content.
- A parameter feeding into a task's script changed — including one you
  consider cosmetic, since the script text is part of the hash.

`SRA_FETCH` uses `storeDir`, so reads specifically are never re-downloaded
even when the task cache misses.

## Sample names look wrong in the VCF

Sample identity comes from the `SM` read-group tag written by `BWA_MEM`, which
comes from the `sample` column of the samplesheet or the `Sample Name` column
of the SRA run table. There is no renaming step to blame. If names are wrong,
fix the input and re-run alignment — check `results/metadata/runs.tsv` to see
what was actually resolved.

## Runs missing from the manifest

`SRA_RESOLVE` reports skipped runs by reason: `not_in_ena`, `not_paired`, or
`organism`. Check the job log. Runs that are single-end, or that ENA lists
without a clean R1/R2 pair, are skipped by design — this pipeline is
paired-end only.
