# Configuration

## Required parameters

| Parameter | Description |
|---|---|
| `--fasta` | Genome assembly FASTA (`.fa`, `.fna`, `.fasta`, or bgzipped `.gz`). Indices, dictionary and intervals are built from it. |
| `--sra_metadata` | NCBI SRA run table. Mutually exclusive with `--input`. |
| `--input` | Samplesheet CSV (`sample,run,fastq_1,fastq_2`). Mutually exclusive with `--sra_metadata`. |

Swapping reference is just a different `--fasta`. There is no genome registry
to maintain, and no assumption that indices already exist.

> **Coordinate compatibility.** Call sets built against different references
> are not comparable and must not be merged. Use `--benchmark_label` to keep
> outputs from different references in clearly separated directories.

## Reads

| Parameter | Default | Description |
|---|---|---|
| `--reads_dir` | `<outdir>/reads` | Where FASTQs live. Existing files are adopted, not re-downloaded. |
| `--download_reads` | `true` | `false` fails the run if anything is missing rather than fetching it. |
| `--max_download_jobs` | `8` | Concurrent downloads. ENA throttles above ~8. |
| `--sra_layout` | `PAIRED` | Layout filter; `ANY` disables it. |
| `--sra_assay` | `null` | e.g. `WGS`. Filters the run table. |
| `--sra_organism` | `null` | e.g. `Cicer reticulatum`. Checked against ENA's `scientific_name`. |

FASTQ naming is fixed by convention: `<RUN>_<SAMPLE>_1.fastq.gz` and
`_2.fastq.gz`. Files already matching that pattern in `--reads_dir` are used
as-is.

## Reference artefacts

| Parameter | Default | Description |
|---|---|---|
| `--fasta_fai` | built | Prebuilt `.fai` |
| `--fasta_dict` | built | Prebuilt `.dict` |
| `--bwa_index` | built | Directory holding a prebuilt bwa index |
| `--chr_regex` | `null` | Only call on sequences whose name matches. `null` means all. |
| `--intervals_min_length` | `0` | Skip sequences shorter than this. |
| `--skip_fasta_normalisation` | `false` | Skip CRLF and header-whitespace cleanup |

### FASTA normalisation

Before anything indexes the assembly, `NORMALISE_FASTA` strips CRLF line
endings and trailing whitespace from header lines. Assemblies arrive in
whatever state they were written in, and both of those quietly break a GATK
stack — a stray CR inside a sequence line shifts coordinates, and header
whitespace makes the `.dict` and the BAM header disagree about sequence names.
`PBA_HatTrick.fasta` ships with CRLF endings, which is what prompted this.

A clean FASTA passes through at no cost. The check reads only the first 2 MB
rather than scanning the whole file.

**Normalising changes byte offsets, which invalidates a `.fai` built against
the original.** Supplying `--fasta_fai`, `--fasta_dict` or `--bwa_index` is
therefore treated as an assertion that the FASTA is already clean, and
normalisation is skipped with a log message rather than silently producing an
index that disagrees with the sequence. For `PBA_HatTrick`, do **not** pass the
`.fai` that ships beside it — let the pipeline normalise and rebuild.

For an assembly with many unplaced scaffolds, `--chr_regex` and
`--intervals_min_length` are the two levers that matter most for cost. The
*C. echinospermum* assembly has 17,304 sequences; scattering across all of them
at 161 samples would create ~2.8 million HaplotypeCaller tasks.

```bash
# PBA_HatTrick: 8 chromosomes + 2 unallocated scaffolds
--chr_regex '_Chr' --intervals_min_length 1000000
```

`PBA_HatTrick` needs neither, strictly — 10 sequences is already a small
scatter — but excluding the two unallocated scaffolds (1.9 Mb of 698.8 Mb)
keeps the call set to placed chromosomes.

## Stage toggles

| Parameter | Default |
|---|---|
| `--trim_reads` | `true` |
| `--skip_qc` | `false` |
| `--skip_markduplicates` | `false` |
| `--skip_variant_calling` | `false` |
| `--skip_filtering` | `false` |
| `--skip_multiqc` | `false` |

## Variant filtering

Hard filters, applied separately to SNPs and indels. VQSR is not used: it
requires a truth set that does not exist for wild *Cicer*.

| Parameter | Default |
|---|---|
| `--snp_filter_expression` | `QD < 2.0 \|\| FS > 60.0 \|\| MQ < 40.0 \|\| MQRankSum < -12.5 \|\| ReadPosRankSum < -8.0 \|\| SOR > 3.0` |
| `--indel_filter_expression` | `QD < 2.0 \|\| FS > 200.0 \|\| ReadPosRankSum < -20.0 \|\| SOR > 10.0` |

These are the GATK recommendations. They are tuned for ~30× human data; at
this cohort's median 7.8× they are a starting point, not a finished answer.
Check the `results/qc/` outputs and the filter summary before trusting them.

## Downstream forks

| Parameter | Default | Description |
|---|---|---|
| `--run_winpca` | `false` | Per-chromosome VCFs for windowed PCA |
| `--run_locus_pca` | `false` | Per-locus VCFs; requires `--loci_bed` |
| `--loci_bed` | `null` | BED of regions; 4th column names each locus |
| `--pca_maf` | `0.05` | Minor allele frequency floor |
| `--pca_max_missing` | `0.10` | Maximum per-site missing fraction |

## Benchmarking

| Parameter | Default | Description |
|---|---|---|
| `--benchmark_label` | timestamp | Tag in trace/report/timeline filenames and the output prefix |
| `--bench_subset` | `null` | Use only the first N samples |
| `--tracedir` | `<outdir>/benchmarks` | Where the evidence lands |

Trace, timeline, report and DAG are always on and never overwrite. See
[Benchmarking](Benchmarking.md).

## Profiles

Combine one execution profile with one software profile.

**Execution**

| Profile | |
|---|---|
| `standard` | Local execution. The default. |
| `uwa` | UWA HPC. SLURM, `work` partition, 96 cores / 1.4 TB ceiling, `/group` and `/scratch` bound into containers. |
| `setonix` | Pawsey Setonix. SLURM, 128 cores / 230 GB per node, work directory forced onto `/scratch`. |
| `slurm` | Generic SLURM with no site assumptions. |

**Software**

`apptainer`, `singularity`, `docker`, `conda`.

**Scaling**

`bench_2`, `bench_8`, `bench_24` set `--bench_subset` for scaling runs.

```bash
-profile uwa,apptainer,bench_8
```

## Resources

Processes are sized by **label**, not individually, so the whole pipeline
re-tunes from `conf/base.config`.

| Label | CPUs | Memory | Time |
|---|--:|--:|--:|
| `process_single` | 1 | 4 GB | 2 h |
| `process_low` | 2 | 8 GB | 4 h |
| `process_medium` | 8 | 32 GB | 8 h |
| `process_high` | 16 | 64 GB | 16 h |
| `process_high_memory` | — | 128 GB | — |
| `process_long` | — | — | 48 h |

Memory scales with `task.attempt`, and tasks killed for OOM (exit 104, 134,
137, 139, 140, 143, 247) retry up to three times with more. At 238 runs with
coverage spanning 0.01× to 146×, a fixed allocation would fail on the
outliers; this is what lets the cohort complete unattended.

Per-process overrides live at the bottom of `conf/base.config`. Override
ceilings with `--max_cpus`, `--max_memory`, `--max_time`.

## Credentials

**Never commit credentials.** `conf/secrets.config` is git-ignored.

The Acacia/AWS keys some setups keep in `~/.nextflow/config` are plaintext.
If that file has been shared, copied, or is on a multi-user system, rotate
the key pair in Acacia.
