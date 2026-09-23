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

## `/usr/bin/env: 'python': No such file or directory` in a GATK task

Exit status 127 from a GATK process under Apptainer:

```
Command error:
  /usr/bin/env: 'python': No such file or directory
```

`gatk` is a Python wrapper script, so it needs `python` on `PATH`.

**Why it happens.** Apptainer inherits the *host* environment by default,
including `PATH` — Docker does not. The `broadinstitute/gatk` image keeps its
interpreter in `/opt/miniconda/envs/gatk/bin`, which is not on the host `PATH`,
so inside the container that directory is never searched and `python` is
invisible. The same image works fine under Docker, which is why this only
shows up on the cluster.

**The fix** is the bioconda build, which installs into `/usr/local/bin` —
a directory that *is* on essentially every host `PATH`, so it resolves
correctly under inherited-environment semantics:

```groovy
container "quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"
```

Same GATK version, different packaging. This is why every container in this
pipeline is a biocontainers image.

The alternative — adding `--cleanenv` to `apptainer.runOptions` so the
container's own `PATH` wins — is not used here: it also strips variables
Nextflow and the scheduler rely on, and trades one class of surprise for
another.

## GenomicsDBImport fails with `IndexOutOfBoundsException: Index: 0`

Look further up the log, not at the stack trace:

```
WARNING  IntervalListCodec  Ignoring interval for unknown reference: <contig>:1-<len>
INFO     IntervalArgumentCollection - Processing 0 bp from intervals
```

Zero intervals resolved, so the column-partition list is empty and
GenomicsDBImport throws on `get(0)`. The exception is a symptom.

**Cause.** A Picard-format `.interval_list` is parsed by `IntervalListCodec`,
which needs an `@SQ` sequence-dictionary header to resolve contig names. A
headerless file does not error — every interval is silently dropped as
"unknown reference".

**Why this is dangerous.** HaplotypeCaller hits the same parse, resolves zero
bases, writes a gVCF containing only headers, and **exits 0**. The task goes
green, the trace records COMPLETED, and the run continues producing nothing
until a downstream step chokes on the emptiness.

**The fix** is plain-text intervals (`.intervals`), one region per line:

```
PBA_HatTrick_Chr2_v1:1-77563191
```

No dictionary requirement, so it cannot fail this way. Measured against the
same reference and contig:

| Format | GATK reports |
|---|---|
| headerless `.interval_list` | `Processing 0 bp from intervals` |
| plain `.intervals` | `Processing 77563191 bp from intervals` |

`BUILD_INTERVALS` emits the plain form and validates every file before
emitting it. `GATK4_HAPLOTYPECALLER` additionally asserts its gVCF contains at
least one record, because a whole chromosome yielding none is never a valid
result. `GATK4_GENOMICSDBIMPORT` now also receives `--reference`, which it
previously lacked entirely.

**If you suspect you have empty gVCFs**, check one directly — exit status will
not tell you:

```bash
zcat sample.chr.g.vcf.gz | grep -vc '^#'     # 0 means no records
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

## Jobs stuck `PENDING (PartitionTimeLimit)`

This is **not** a queue wait — it is permanent. The job asks for more wall time
than the partition allows, so SLURM will never schedule it:

```
     JOBID   NAME              PARTITION  TIME_LIMIT  STATE    REASON
       785   nf-ALIGN_READS_B  benchmarki 3-00:00:00  PENDING  (PartitionTimeLimit)
```

Compare what the job asks for against the partition maximum:

```bash
sinfo -o "%.14P %.12l %.6D %.20C"        # TIMELIMIT column
squeue -u $USER -o "%.10i %.16j %.10P %.11l %R"
```

On UWA:

| Partition | Nodes | Wall limit |
|---|--:|---|
| `work` | 12 | 3 days |
| `benchmarking` | 1 (k003) | **24 h** |
| `ondemand` | 1 | 12 h |

`conf/base.config` requests 72 h for `BWA_MEM`, `GENOMICSDBIMPORT` and
`GENOTYPEGVCFS`, which is fine on `work` and impossible on `benchmarking`.

**The fix** is the `benchmarking` profile, which lowers `params.max_time` to
24 h so `check_max()` caps every request to fit:

```bash
sbatch bin/run_pipeline.sbatch -profile uwa,benchmarking,apptainer,bench_8 ...
```

Do not simply raise `max_time` back up to make the numbers look right — the
partition limit is real, and the job will go straight back to
`PartitionTimeLimit`.

### What fits on the benchmarking node, and what does not

The node is uncontended, which is exactly what you want for measurements. But
24 h is a real constraint:

- **Alignment and HaplotypeCaller fit.** At ~8x median on a ~700 Mb genome,
  `BWA_MEM` is hours.
- **The deepest samples may not.** `Besev_079` is ~146x, about 19x the median.
  It is not in `bench_2` or `bench_8`.
- **Joint genotyping at full cohort size will not.** `GENOTYPEGVCFS` also wants
  512 GB, a third of the node. Run that stage on `work`.

### Running everything on `benchmarking`

Running the whole pipeline there is a reasonable first move — the node is
uncontended, and 512 GB fits comfortably in its 1.5 TB. Only the 24 h ceiling
is in question, and the honest way to find out whether it binds is to run and
see.

The profile is set up so a timeout tells you clearly rather than quietly
burning days:

```groovy
// benchmarking.config
errorStrategy = { task.exitStatus in [104,134,137,139,247] ? 'retry' : 'finish' }
```

Out-of-memory (137 and friends) still retries, because memory doubles per
attempt and the node has headroom. **A timeout does not retry.** `max_time` is
already pinned at the partition ceiling, so every retry would get the same
24 h and fail identically — four attempts would spend four days learning
nothing.

If a stage does time out, the fix is a partition with a longer limit, not
another attempt:

```bash
# that stage only, on work (3-day limit)
sbatch bin/run_pipeline.sbatch -profile uwa,apptainer --outdir results/full -resume
```

`-resume` keeps everything already computed on the benchmarking node.

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

## `~/.nextflow/config` silently affects every run

Nextflow merges `~/.nextflow/config` into **every** run, before the project's
own config. Two consequences worth knowing about.

**Credentials leak into runs you did not intend.** A top-level `aws { accessKey
= '...' secretKey = '...' }` block there is loaded by every pipeline you launch
on that account, not just the one it was added for. Check with:

```bash
nextflow config -profile uwa | grep -A6 '^aws'
```

Keep credentials in a file you pass explicitly instead, and rotate anything
that has been sitting in a plaintext home config:

```bash
cp conf/secrets.config.template conf/secrets.config   # git-ignored
chmod 600 conf/secrets.config
sbatch bin/run_pipeline.sbatch -profile uwa,apptainer -c conf/secrets.config ...
```

**Profile names collide.** If your home config defines a profile called
`slurm`, `conda` or `apptainer`, it merges with this project's profile of the
same name. That is usually harmless and occasionally not — a home `slurm`
profile setting `process.queue = "work"` will apply on a cluster with no such
partition.

To see what a run will actually use:

```bash
nextflow config -profile uwa,apptainer          # merged: home + project
nextflow config -C nextflow.config -profile uwa # project only
```

`-C` (capital) ignores every other config source. It is useful for diagnosis
but not for real runs, since it also discards `-params-file` and site settings.

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
