# Pipeline steps

## 1. Resolve the run table — `SRA_RESOLVE`

Reads an NCBI SRA run table and resolves each run to its paired FASTQ URLs,
MD5 checksums and byte counts via the ENA portal API.

`bin/parse_runtable.py` handles what NCBI actually exports: UTF-16 or UTF-8,
comma or tab separated, with or without a BOM, and with the several column
names NCBI uses for the same field across export paths.

**On the ENA query.** The portal's `filereport` endpoint returns an *empty
result set*, not an error, for a comma-separated list of accessions. So the
parser resolves by study instead: it reads the `BioProject`/`SRA Study` column
and fetches that study once, which covers any cohort size in a single request.
Where no study is named it falls back to the `search` endpoint with batched
`run_accession="X" OR ...` queries.

Filters (`--sra_layout`, `--sra_assay`, `--sra_organism`) are applied before
the lookup. Runs that resolve to something other than a clean R1/R2 pair are
reported and skipped rather than silently dropped — ENA sometimes lists an
unpaired "orphan" file alongside the mates, so the parser selects by filename
suffix rather than by position in the list.

**Output:** `results/metadata/runs.tsv`

## 2. Fetch reads — `SRA_FETCH`

Downloads each pair and verifies both MD5s.

`storeDir` points at `--reads_dir`, so a run whose files are already there
never executes. That makes re-running free, and lets a download started with
`bin/fetch_reads.sh` be picked up mid-flight.

A file failing checksum is deleted and the task fails, so a corrupt transfer
is retried rather than silently poisoning the alignment. `maxForks` is bound
by `--max_download_jobs` because ENA throttles aggressively — measured
throughput at 16 concurrent connections was *half* that at 8.

## 3. Prepare the reference — `PREPARE_GENOME`

From the FASTA alone, builds and publishes to `results/reference/`:

| Step | Output |
|---|---|
| `SAMTOOLS_FAIDX` | `.fai` (and `.gzi` for bgzipped input) |
| `GATK4_CREATESEQUENCEDICTIONARY` | `.dict` |
| `BWA_INDEX` | `bwa/` |
| `BUILD_INTERVALS` | one `.interval_list` per sequence |

Each artefact is resolved **by name**. The original destructured a glob
positionally:

```groovy
map { bed, fai, pac, sa, amb, genome, ann, bwt, fa, dict -> ... }
```

which binds each index to whatever the glob happened to sort into that slot.
Adding or removing a single sidecar file silently mis-assigns every reference
input, with no error.

`BUILD_INTERVALS` emits one interval list per sequence, filtered by
`--chr_regex` and `--intervals_min_length`, and fails loudly if nothing
matches. The original derived interval names with `splitText()`, which keeps
the trailing newline on each line and passed it straight into GATK's `-L`.

## 4. Read QC and trimming — `FASTQC`, `FASTP`

fastp does adapter detection and quality trimming; FastQC reports on the raw
reads. Both feed MultiQC. Disable with `--skip_qc` / `--trim_reads false`.

## 5. Align — `BWA_MEM`

`bwa mem` piped into `samtools sort`. One thread is reserved for the sort so
alignment and sorting overlap; sort memory is derived from `task.memory`.

**The read group is the important part:**

```
@RG  ID:<run>  SM:<sample>  PL:ILLUMINA  LB:<sample>  PU:<run>
```

`ID` is the sequencing run, `SM` is the biological sample. This single
decision is what removes the original's `PullSampleName` → `RenameVCFs` pair:
sample identity is correct in the BAM from the first step, so nothing
downstream has to repair it.

## 6. Merging runs into samples — `SAMTOOLS_MERGE`

**This is required for this cohort, not an optimisation.** The 238 runs are
161 biological samples:

| Runs per sample | Samples |
|--:|--:|
| 1 | 85 |
| 2 | 75 |
| 3 | 1 |

After alignment the workflow re-keys each BAM by `meta.sample`, groups, and
branches. Single-run samples pass straight through; multi-run samples are
merged. Because the RG `SM` tag already matches across a sample's runs, the
merged BAM is a valid single-sample BAM with no header surgery.

Treating runs as samples would give you 238 "samples", 77 of which are partial
duplicates of another — which would badly distort any PCA.

## 7. Mark duplicates — `GATK4_MARKDUPLICATES`

Java heap is derived from `task.memory` (80%, leaving headroom for
off-heap and native allocation). The original hardcoded `-Xmx4g` inside a
64 GB request, which both wasted 60 GB of the allocation and forced Picard's
sorting collection to spill to disk unnecessarily.

Produces `<sample>.bam` and `<sample>.bam.bai` in `results/alignments/`.

## 8. Alignment QC — `SAMTOOLS_STATS`

`samtools stats`, `flagstat` and `coverage` per sample. Beyond QC, realised
coverage is what converts a wall-clock measurement into a per-Gbp cost model
for the allocation case.

## 9. Joint genotyping — `CALL_VARIANTS`

```
HAPLOTYPECALLER    sample × interval          →  gVCF shard
GENOMICSDBIMPORT   interval, all samples      →  GenomicsDB workspace
GENOTYPEGVCFS      interval                   →  cohort VCF shard
BCFTOOLS_CONCAT    all shards                 →  genome-wide cohort VCF
```

**Why this order.** The original combined gVCFs per *sample* across intervals
before importing. That serialises work the scatter had just parallelised, and
produces an intermediate the cohort join does not need. Here nothing crosses
the cohort except the two per-interval joins.

`GENOMICSDBIMPORT` is the pipeline's largest memory consumer and usually the
step that dictates an allocation request. It builds a sample map rather than a
command line of `-V` arguments, so cohort size does not run into argument
limits, and `--batch-size 50` bounds open file handles — with 161 samples you
exhaust file descriptors well before memory.

`BCFTOOLS_CONCAT` orders shards by their position in the `.fai`, so the output
is coordinate-sorted without a separate sort.

## 10. Filter — `GATK4_VARIANTFILTRATION`

SNPs and indels are separated, filtered with their own expressions, and merged
back. See [Configuration](Configuration.md#variant-filtering) — the defaults
are GATK's human-oriented recommendations and warrant checking against this
cohort's shallow coverage.

## 11. Downstream forks — `DOWNSTREAM`

See [Downstream analysis](Downstream-Analysis.md).
