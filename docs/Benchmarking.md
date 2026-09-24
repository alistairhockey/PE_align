# Benchmarking

This pipeline exists partly to produce the evidence for a Pawsey allocation
application. Instrumentation is always on and never overwrites.

## What is captured

Every run writes four files to `<outdir>/benchmarks/`, tagged with
`--benchmark_label` (or a timestamp):

| File | Contents |
|---|---|
| `trace-<tag>.txt` | One row per task: CPU%, peak RSS, peak VMEM, wall time, I/O, exit status, retries |
| `report-<tag>.html` | Interactive summary with per-process distributions |
| `timeline-<tag>.html` | Gantt chart — shows where parallelism stalls |
| `dag-<tag>.html` | Workflow graph |

The trace carries the full field set, including `peak_rss`, `peak_vmem`,
`rchar`, `wchar`, `read_bytes`, `write_bytes` and `attempt`. Those are what an
allocation case is actually built from; the defaults omit most of them.

## The scaling protocol

**Do not extrapolate from a single run, and do not start at 238.** Run three
increasing subsets, then project.

```bash
for n in 2 8 24; do
  sbatch bin/run_pipeline.sbatch -profile uwa,apptainer,bench_${n} \
      --sra_metadata assets/cret_metadata.tsv \
      --fasta /group/peg/cicer/chickpea/genome/PBA_HatTrick/PBA_HatTrick.fasta \
      --reads_dir /group/peg/cicer/cret/reads \
      --outdir results_bench${n} \
      --benchmark_label bench${n}
done
```

Subsetting is applied **per sample, never per run**, so a sample never loses
half of its data and its coverage stays honest.

Three points let you separate the two scaling regimes:

- **Per-sample stages** (`FASTP`, `BWA_MEM`, `MARKDUPLICATES`,
  `HAPLOTYPECALLER`) scale linearly in sample count. Extrapolation is safe.
- **Cohort-wide stages** (`GENOMICSDBIMPORT`, `GENOTYPEGVCFS`) scale
  super-linearly in memory, because every sample's gVCF for an interval is
  held at once. These must be *measured* at increasing cohort size, not
  extrapolated from one point. This is the single biggest risk in the
  application: a memory figure derived from 8 samples will understate 161.

## Analysing a run

```bash
bin/summarise_benchmark.py results/benchmarks/trace-bench8.txt
```

```
process                             n  cpu   cpu%  med time  max time   req mem   p95 rss   mem%         rec
BWA_MEM                            20   16    77%     1.42h     2.47h   48.0 GB   15.5 GB    32%     16c/20G
GATK4_GENOMICSDBIMPORT             12    4    23%     2.50h     5.49h   64.0 GB   40.5 GB    63%      1c/51G
GATK4_HAPLOTYPECALLER              30    4    56%     1.03h     1.95h   16.0 GB    8.6 GB    54%      4c/11G
```

Read it as:

- **`cpu%`** — mean utilisation against *allocated* cores. `BWA_MEM` at 77% is
  healthy. `GENOMICSDBIMPORT` at 23% means three of its four cores idle, and
  the recommendation drops it to 1.
- **`mem%`** — p95 peak RSS against the request. `BWA_MEM` at 32% is
  over-requesting by 3×.
- **`rec`** — right-sized `cpus/memory` from the p95 plus 25% headroom.

Project to the full cohort:

```bash
bin/summarise_benchmark.py results/benchmarks/trace-bench8.txt \
    --scale-from 8 --scale-to 161 --format markdown -o benchmark.md
```

## Service units

```bash
bin/estimate_su.py results/benchmarks/trace-bench8.txt \
    --scale-from 8 --scale-to 161
```

Setonix CPU nodes are dual AMD EPYC 7763: **128 cores, ~230 GB usable**, so
about **1.79 GB per core**. That ratio drives the whole cost model.

**Memory-implied cores.** A task requesting 64 GB occupies roughly 36 cores'
worth of a node whether or not it uses them, because no one else can use that
memory. The estimator reports this as `mem-cores` and marks the task's binding
as `memory` or `cores`:

```
process                  cores   req mem  mem-cores    bind      core-h      SU core
GATK4_GENOMICSDBIMPORT       4       64G       35.6  memory      27,104       27,104
BWA_MEM                     16       48G       26.7  memory      18,104       18,104
```

Both are memory-bound. Their core counts are nearly irrelevant to what they
cost. **Right-sizing memory on memory-bound processes is the cheapest way to
reduce the allocation** — and, because it comes from measurement, the easiest
to defend.

Two accounting models are reported. Whole-node is modelled as core-hours
divided by a packing efficiency (default 80%), *not* a whole node per task —
the scheduler packs many tasks onto each node, and charging per task
overstates the total several-fold.

> Charge rate, node geometry and packing efficiency are all flags, not
> hardcoded policy. **Confirm the current rate and partition accounting with
> Pawsey before quoting any figure in an application:**
> `--su-per-core-hour`, `--node-cores`, `--node-mem-gb`,
> `--packing-efficiency`.

An `--overhead` multiplier (default 1.15) covers retries, requeues and the
reruns any real project needs.

## Where the resource figures come from

The attempt-1 values in `conf/base.config` are **measured, not estimated**.
They derive from settings tested against the *C. echinospermum* cohort on a
comparable ~700 Mb *Cicer* genome, where they were found sufficient to avoid
OOM:

| Process | Attempt 1 | Time |
|---|--:|--:|
| `BWA_MEM` | 16 cpus / 128 GB | 72 h |
| `GATK4_MARKDUPLICATES` | 4 cpus / 64 GB | 24 h |
| `SAMTOOLS_MERGE` | 4 cpus / 64 GB | 24 h |
| `GATK4_HAPLOTYPECALLER` | 4 cpus / 64 GB | 24 h |
| `GATK4_GENOMICSDBIMPORT` | 4 cpus / 64 GB | 72 h |
| `GATK4_GENOTYPEGVCFS` | 4 cpus / **512 GB** | 72 h |

Memory doubles per retry on top of these, so they are a floor rather than a
ceiling. Do not reduce them without a trace to justify it — `p95 peak RSS` from
`summarise_benchmark.py` is the right basis for tightening, and tightening is
what makes an allocation request defensible. Guessing lower is what produces
OOM.

### This dictates the Setonix partition

`GATK4_GENOTYPEGVCFS` at 512 GB **does not fit on a standard Setonix node**:

| Partition | Cores | Memory |
|---|--:|--:|
| `work` / `long` | 128 | ~230 GB |
| `highmem` | 128 | ~980 GB |

A 512 GB request capped to 230 GB would simply OOM again, so `conf/setonix.config`
sets `max_memory` to the highmem node and routes by what each task needs:

```groovy
queue = {
    task.memory && task.memory > 230.GB ? 'highmem'
  : task.time   && task.time   > 24.h   ? 'long'
  :                                       'work'
}
```

**This is a load-bearing point for the allocation application.** Joint
genotyping this cohort requires `highmem`, which is charged at a higher rate
than `work`. Budget for it explicitly rather than discovering it mid-run, and
confirm the current highmem charge rate with Pawsey — `bin/estimate_su.py`
takes `--node-mem-gb` and `--su-per-core-hour` as flags for exactly this.

## The allocation application table

`bin/allocation_table.py` groups processes into the stages a reviewer thinks in
and keeps the arithmetic self-consistent — **CPU h = Jobs x Cores/job x
Wall-time** — because that is the first thing anyone checks.

```bash
bin/allocation_table.py --since 2026-09-20                              # measured
bin/allocation_table.py --since 2026-09-20 --scale-from 24 --scale-to 161
bin/allocation_table.py --since 2026-09-20 --format markdown -o table.md
```

Column meanings, which matter if you are asked to defend them:

| Column | What it is |
|---|---|
| Jobs | Task count |
| Cores/job | Cores weighted by wall time, so a stage mixing 16-core alignment with 1-core indexing reports the cores that actually cost something |
| Wall-time (h) | **Mean per job**, so the CPU h product holds |
| Memory (GB) | Observed peak RSS + 25% headroom — what is worth requesting, not what was requested |
| CPU h | Jobs x Cores/job x Wall-time |

Stages with no tasks print an em dash rather than a zero. Do not fill those in
by hand: an unrun stage has no measurement, and a reviewer is entitled to ask
where a number came from.

## Coverage heterogeneity

This cohort is not uniform, and a mean is misleading:

| | Coverage (740 Mb genome) |
|---|---|
| Median | 7.8× |
| IQR | 5.3× – 12.3× |
| Range | 0.01× – 145.8× |
| Below 5× | 30 samples |
| Below 8× | 82 samples |

Consequences for the benchmark:

1. **Pick subsets that span the range.** A `bench_8` of only shallow samples
   will understate `BWA_MEM` badly. `Besev_079` alone is 145× — nearly 19× the
   median — and will dominate any subset it appears in.
2. **`CudiA_122` has 9.5 Mbp total (~0.01×).** It will produce a BAM and a
   gVCF of almost entirely no-calls. Decide whether to exclude it before the
   full run rather than discovering it in the PCA.
3. **Shallow coverage inflates genotype missingness**, which is the dominant
   artefact in downstream PCA. See
   [Downstream analysis](Downstream-Analysis.md).

## Right-sizing, then re-measuring

The loop is: measure → apply the recommendations to `conf/base.config` →
re-run a subset → confirm nothing now OOMs. A request derived this way is both
cheaper and far more credible than a padded one.

Keep memory headroom above the p95 rather than the median: the retry ladder in
`conf/base.config` exists precisely because outliers are real, but a retry
costs a full task's wall time and you do not want it firing routinely.
