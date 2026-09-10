# Cr_align

Alignment and joint variant calling for *Cicer reticulatum* whole-genome
short-read data, instrumented for resource benchmarking.

Give it an NCBI SRA run table and a genome assembly. It resolves the runs
against ENA, downloads the reads, builds every reference index it needs,
aligns, marks duplicates, joint-genotypes the cohort, and writes a full
per-task record of CPU and memory use.

```bash
sbatch bin/run_pipeline.sbatch -profile uwa,apptainer \
    --sra_metadata assets/cret_metadata.tsv \
    --fasta /path/to/assembly.fna \
    --reads_dir /group/peg/cicer/cret/reads \
    --outdir results
```

Nothing runs on a login node: the sbatch wrapper puts the Nextflow driver in a
small allocation, and the driver submits every task as its own job.

Everything else — FASTQ downloads, `.fai`, `.dict`, the bwa index, scatter
intervals, the samplesheet — is derived.

---

## What it produces

| Output | Path |
|---|---|
| Per-sample BAM + index | `results/alignments/` |
| Cohort VCF, raw and filtered | `results/variants/` |
| Per-interval VCF shards | `results/variants/per_interval/` |
| Reference indices, intervals | `results/reference/` |
| QC (FastQC, fastp, samtools, MultiQC) | `results/qc/` |
| PCA-ready SNP matrix | `results/pca/input/` |
| Per-chromosome VCFs (winPCA) | `results/pca/winpca/` |
| Per-locus VCFs | `results/pca/loci/` |
| Trace, timeline, report, DAG | `results/benchmarks/` |

## The cohort

238 paired-end WGS runs from BioProject **PRJNA416007** (SRP123332),
Illumina HiSeq 4000, all *Cicer reticulatum*.

Those 238 runs are **161 biological samples**: 85 samples have one run,
75 have two, and 1 has three. Per-sample BAM merging is therefore required,
not optional — the pipeline keeps `sample` and `run` distinct in read-group
metadata and merges on `sample`.

Coverage is heterogeneous and mostly shallow:

| | Coverage (assuming a 740 Mb genome) |
|---|---|
| Median | 7.8× |
| Interquartile range | 5.3× – 12.3× |
| Range | 0.01× – 145.8× |
| Samples below 5× | 30 |
| Samples below 8× | 82 |

`CudiA_122` has 9.5 Mbp total (~0.01×) and cannot be genotyped meaningfully;
`Besev_079` has 145× and dominates alignment cost. Both matter for
interpreting the benchmark and for deciding a minimum-coverage cutoff before
PCA. See the [wiki](docs/) for the full discussion.

## Requirements

- Nextflow ≥ 24.04
- One of: Apptainer/Singularity, Docker, or Conda
- ~0.71 TB for the FASTQs, plus working space for BAMs and gVCFs

## Documentation

The `docs/` directory is the wiki source.

| Page | |
|---|---|
| [Quick start](docs/Quick-Start.md) | First run, smallest possible |
| [Configuration](docs/Configuration.md) | Every parameter and profile |
| [Pipeline steps](docs/Pipeline-Steps.md) | What each process does and why |
| [Benchmarking](docs/Benchmarking.md) | Producing the allocation case |
| [Downstream analysis](docs/Downstream-Analysis.md) | winPCA and locus PCA forks |
| [Troubleshooting](docs/Troubleshooting.md) | Failure modes and fixes |
| [Changes from the original](docs/Changes-From-Original.md) | What was fixed and why |

## Layout

```
main.nf                     workflow entry point and parameter validation
nextflow.config             parameters, profiles, benchmarking instrumentation
conf/
  base.config               resources by process label
  uwa.config                UWA HPC (SLURM)
  setonix.config            Pawsey Setonix
  test.config               minimal end-to-end test
modules/local/              one process per file
subworkflows/local/         fetch, prepare genome, align, call, downstream
bin/                        helper scripts (see below)
assets/                     cohort metadata and resolved run manifest
docs/                       wiki source
```

| Script | Purpose |
|---|---|
| `bin/parse_runtable.py` | SRA run table → resolved manifest with URLs and MD5s |
| `bin/fetch_run.py` | Download one FASTQ and verify its checksum |
| `bin/fetch_reads.sh` | Bulk download outside Nextflow, resumable |
| `bin/fetch_reads_chain.sh` | Submit that download as backfill-friendly chained jobs |
| `bin/run_pipeline.sbatch` | Submit the Nextflow driver to SLURM |
| `bin/stop_fetch.sh` | Stop a bulk download cleanly |
| `bin/summarise_benchmark.py` | Trace → per-process resource table and right-sizing |
| `bin/estimate_su.py` | Resource table → Setonix service-unit projection |

## Licence

MIT. See [LICENSE](LICENSE).
