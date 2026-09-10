# Changes from the original

The original pipeline is preserved at commit `5aff337` ("baseline — original
pipeline as inherited"), so every change here is diffable against its starting
point:

```bash
git diff 5aff337 HEAD -- main.nf
```

## Correctness

### `groupTuple(by: 3)` indexed past the end of the tuple

```groovy
| combine(chromosomeNames)
| map { id, vcf, chrName -> [ id, vcf, clean_chrName ] }
| groupTuple(by: 3) | view
```

The tuple has three elements — indices 0, 1, 2. Grouping by index 3 is out of
range. The pipeline ended on `| view` here, with `ConsolidateGVCFs` and
`JointCallCohort` commented out, so **it never produced a joint-called VCF**.
The scatter-gather rewrite removes this path entirely.

### The reference channel was bound by glob sort order

```groovy
Channel.fromPath("${params.genome}.*")
| collect
| map { bed, fai, pac, sa, amb, genome, ann, bwt, fa, dict -> ... }
```

This destructures a *glob result* positionally, so each index file is bound to
whatever alphabetical sorting happened to put in that slot. Adding or removing
a single sidecar file — a `.bed`, a stray `.gzi` — shifts every subsequent
binding, silently handing bwa the wrong file with no error.

`PREPARE_GENOME` now resolves each artefact **by name** and builds anything
missing.

### `HaplotypeCaller` ended on a dangling line continuation

```groovy
  -L ${intName} \\  
```

A trailing backslash followed by whitespace and end-of-script. Combined with
interval names that still carried their newline (below), the emitted command
was fragile at best.

### `splitText()` kept trailing newlines

```groovy
Channel.fromPath(params.intervals) | splitText()
```

`splitText()` retains the line terminator, so each interval name reached GATK
as `"chr1\n"`. `BUILD_INTERVALS` now emits proper `.interval_list` files, one
per sequence, derived from the `.fai`.

### Sample renaming worked around a problem that should not exist

The original set `SM` to a value that later needed correcting, so it ran
`bcftools query -l` → `picard RenameSampleInVcf` on every gVCF
(`PullSampleName` → `RenameVCFs`).

`BWA_MEM` now writes `ID:<run> SM:<sample>` at alignment. Sample identity is
correct from the first step, and both processes are deleted.

## Efficiency

### Per-sample `CombineGVCFs` serialised the scatter

The original scattered `HaplotypeCaller` by interval, then immediately
combined each sample's shards back together *before* the cohort join —
undoing the parallelism it had just created, to produce an intermediate
`GenomicsDBImport` does not need.

Now: `HaplotypeCaller` per sample × interval → `GenomicsDBImport` per interval
across all samples → `GenotypeGVCFs` per interval → `bcftools concat`. Nothing
crosses the whole cohort except the two per-interval joins.

### Java heap did not match the allocation

```groovy
process MarkDuplicates {
  memory '64GB'
  """
  picard -Xmx4g MarkDuplicates ...
  """
}
```

64 GB requested, 4 GB usable. The other 60 GB was reserved and idle while
Picard's sorting collection spilled to disk. Heap is now derived from
`task.memory` in every Java process.

### Resources were declared ad hoc

Several processes declared `memory '64 GB'` with no `cpus`; `Indexation`,
`IdxMerge` and others declared nothing at all and silently took the default.
Resources are now assigned by **label** in `conf/base.config`, with
per-process overrides in one place, so the pipeline re-tunes from a single
file after benchmarking.

### `MergeBAMs` wrote an intermediate it did not declare

```groovy
samtools merge -o ${accession}.merged.bam ${bams}
samtools sort ${accession}.merged.bam > ${accession}.bam
```

Both files stayed in the work directory. `bwa mem` already emits
coordinate-sorted BAMs, so the merge output needs no re-sort.

## Portability

### Every path was hardcoded to a filesystem that no longer exists

```groovy
genome       = '/scratch/sae001/cicer-data/cechi/genome/...'
modules_path = '/scratch/sae001/ahockey/all_ce_align/modules.nf'
publishDir '/scratch/sae001/ahockey/all_ce_align/alignments'
```

`/scratch/sae001` is not present on this machine. `modules_path` as a
*parameter* also meant the pipeline could not be moved without editing it, and
`include ... from params.modules_path` defeats the module system.

Now: standard relative includes, `publishDir` derived from `--outdir`, and
every path a parameter.

### Two inputs, everything else derived

The original required a pre-built reference index set, pre-staged reads and
a matching interval file. It now takes an SRA run table and a FASTA, and
builds the rest — including downloading and checksum-verifying the reads.

Reference selection is just a different `--fasta`; there is no genome registry
to maintain.

### No retry, no error strategy

A single OOM killed the run. With coverage spanning 0.01× to 145.8×, that is
not hypothetical. Memory now scales with `task.attempt` and OOM/kill exit
codes retry up to three times.

### Containers were unpinned and half-specified

```groovy
container = 'broadinstitute/gatk'   // no tag: whatever :latest is today
```

Some processes had a `conda` directive and no container, some both, some
neither. Every process now pins an exact version for both, and bioconda images
carry the `quay.io/` registry they actually live on.

## Additions

| Area | |
|---|---|
| Benchmarking | Trace, timeline, report and DAG always on, with the full resource field set; `summarise_benchmark.py` and `estimate_su.py` |
| QC | FastQC, fastp, samtools stats/flagstat/coverage, MultiQC |
| Filtering | GATK hard filters applied separately to SNPs and indels |
| Downstream | PCA-ready SNP matrix, per-chromosome and per-locus forks |
| Validation | Parameters checked up front with actionable messages; duplicate sample+run pairs and unpaired rows rejected |
| Scheduler | `sbatch` wrappers; no step runs on a login node |
| Reproducibility | Every process emits a `versions.yml` |

## Things deliberately not carried over

- **The single-end branch.** The original mixed paired and single-end reads in
  one channel. This cohort is entirely paired-end WGS, and the mixed handling
  was a source of confusion; unpaired rows are now rejected explicitly.
- **`ConsolidateGVCFs` / `JointCallCohort` as written.** Replaced by
  `GATK4_GENOMICSDBIMPORT` and `GATK4_GENOTYPEGVCFS` with a sample map,
  bounded batch size and per-interval scatter.
- **The Acacia/Fusion profile from `~/.nextflow/config`.** It contains
  plaintext credentials and must not enter the repository. Keep site
  credentials in `conf/secrets.config`, which is git-ignored. **If those keys
  have been shared or sat on a multi-user system, rotate them in Acacia.**
