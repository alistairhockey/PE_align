# PE_align

Alignment and joint variant calling for *Cicer reticulatum* whole-genome
short-read data, with per-task resource accounting for HPC allocation
applications.

The pipeline takes **two inputs** — an NCBI SRA run table and a genome
assembly — and derives everything else.

```bash
nextflow run . -profile uwa,apptainer \
    --sra_metadata assets/cret_metadata.tsv \
    --fasta /group/peg/cicer/chickpea/genome/PBA_HatTrick/PBA_HatTrick.fasta \
    --reads_dir /group/peg/cicer/cret/reads \
    --outdir results
```

## Pages

| Page | What it covers |
|---|---|
| [Quick start](Quick-Start.md) | Smallest working run, then the real one |
| [Configuration](Configuration.md) | Every parameter, profile and resource label |
| [Pipeline steps](Pipeline-Steps.md) | Each process, what it does and why |
| [Benchmarking](Benchmarking.md) | Scaling runs, right-sizing, SU projection |
| [Downstream analysis](Downstream-Analysis.md) | winPCA and locus-specific PCA forks |
| [Troubleshooting](Troubleshooting.md) | Failure modes and how to recover |
| [Changes from the original](Changes-From-Original.md) | Every fix, and the bug it addressed |

## Design in one page

```
SRA run table ──► SRA_RESOLVE ──► manifest (URL + MD5 per file)
                                       │
                                       ▼
                                  SRA_FETCH  ◄── reuses --reads_dir
                                       │
genome FASTA ──► PREPARE_GENOME        │
                  ├─ SAMTOOLS_FAIDX    │
                  ├─ CREATESEQUENCEDICT│
                  ├─ BWA_INDEX         │
                  └─ BUILD_INTERVALS   │
                        │              │
                        ▼              ▼
                    ALIGN_READS: FASTP ─► BWA_MEM ─► merge by sample
                                              ─► MARKDUPLICATES ─► BAM
                                                        │
                              ┌─────────────────────────┘
                              ▼
   CALL_VARIANTS:  HAPLOTYPECALLER   (sample × interval)
                   GENOMICSDBIMPORT  (interval, all samples)
                   GENOTYPEGVCFS     (interval)
                   BCFTOOLS_CONCAT   ─► cohort VCF
                   VARIANTFILTRATION ─► filtered cohort VCF
                              │
                              ▼
   DOWNSTREAM:     VCF_PREP_PCA ─┬─► VCF_SPLIT_CHROM  (--run_winpca)
                                 └─► VCF_SUBSET_LOCI  (--run_locus_pca)
```

Two principles shape the structure:

**Sample identity is set once, at alignment.** The read group written by
`BWA_MEM` carries `ID` = run and `SM` = biological sample. Everything
downstream reads sample identity from the BAM, so there is no renaming step
and no opportunity for names to drift.

**Scatter early, join once.** `HaplotypeCaller` runs per sample per interval.
The only cohort-wide operations are `GenomicsDBImport` and `GenotypeGVCFs`,
both of which run per interval. Nothing serialises across the whole cohort.

## The data

238 paired-end WGS runs, BioProject `PRJNA416007`, all *Cicer reticulatum*,
Illumina HiSeq 4000 — but only **161 biological samples**. See
[Pipeline steps](Pipeline-Steps.md#merging-runs-into-samples) for why that
distinction drives the whole design, and
[Benchmarking](Benchmarking.md#coverage-heterogeneity) for what the coverage
spread means for resource planning.
