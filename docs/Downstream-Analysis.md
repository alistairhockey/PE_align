# Downstream analysis

The pipeline stops at **analysis-ready inputs**, not at finished analyses.
Windowed and locus-specific PCA scripts stay external so they can be iterated
on without re-running variant calling — which, at this cohort size, is the
difference between a minute and a week.

## The fork points

```
filtered cohort VCF
        │
        ▼
   VCF_PREP_PCA                    results/pca/input/
   biallelic PASS SNPs,
   MAF and missingness filtered
        │
        ├──► VCF_SPLIT_CHROM       results/pca/winpca/   (--run_winpca)
        │    one VCF per chromosome
        │
        └──► VCF_SUBSET_LOCI       results/pca/loci/     (--run_locus_pca)
             one VCF per named region
```

## Preparing the SNP matrix — `VCF_PREP_PCA`

```bash
--run_winpca --pca_maf 0.05 --pca_max_missing 0.10
```

Keeps biallelic SNPs with a `PASS` filter, then applies a minor-allele
frequency floor and a per-site missingness ceiling.

**Why this matters more than usual here.** This cohort has a median coverage
of 7.8×, and 82 of 161 samples are below 8×. At that depth, genotype
missingness correlates strongly with sequencing depth — so if you do not
control it, **PC1 will track per-sample depth rather than population
structure**. Shallow samples cluster with each other because they share
missing sites, not ancestry.

Practical consequences:

- `--pca_max_missing 0.10` is a reasonable default, but check
  `results/pca/input/*.stats.txt` for how many sites survive. If the count
  collapses, your cohort cannot support that threshold and you should either
  relax it or drop the shallowest samples.
- Consider excluding samples below ~3–5× entirely. `CudiA_122` (0.01×) should
  certainly go.
- Plot PC1 against per-sample mean depth from `results/qc/samtools/*.coverage`
  before interpreting anything. If they correlate, tighten the filters or drop
  samples — do not interpret the axes.

Written to `results/pca/input/<prefix>.pca_input.vcf.gz` with a stats file
recording sites in, sites out, and sample count.

## Windowed PCA — `--run_winpca`

Produces one VCF per chromosome in `results/pca/winpca/`, named
`<prefix>.<chromosome>.vcf.gz`, plus `chrom_list.txt`.

Per-chromosome inputs let a windowed PCA parallelise across chromosomes and
keep memory bounded — a genome-wide VCF for 161 samples is awkward to slide a
window over in one process.

Restrict to real chromosomes rather than scaffolds:

```bash
--chr_regex '^cicec\.S2Drd065\.gnm1\.chr'
```

Without it, an assembly with 17,304 sequences produces 17,304 VCFs.

## Locus-specific PCA — `--run_locus_pca`

```bash
--run_locus_pca --loci_bed regions.bed
```

`regions.bed` is a standard BED whose **4th column names each locus**:

```
cicec.S2Drd065.gnm1.chr3	12400000	14750000	INV1
cicec.S2Drd065.gnm1.chr1	3100000	3450000	FT_region
```

Without that name column the outputs would be labelled by coordinates and the
resulting plots would be unreadable. BED is 0-based half-open; the module
converts to bcftools' 1-based inclusive regions for you.

Produces `results/pca/loci/<prefix>.<locus>.vcf.gz` and `loci_summary.tsv`
recording sites and sample count per locus — check that summary before
plotting, since a locus with few surviving sites will give a meaningless PCA.

## Hooking up your own scripts

Point them at the fork outputs:

| Analysis | Input |
|---|---|
| Windowed PCA | `results/pca/winpca/*.vcf.gz` |
| Locus PCA | `results/pca/loci/*.vcf.gz` |
| Anything genome-wide | `results/pca/input/*.pca_input.vcf.gz` |

When you are ready to bring those scripts into the pipeline, the pattern is:
add a module under `modules/local/`, include it in
`subworkflows/local/downstream.nf`, and gate it behind a `--run_*` flag.
Keeping each analysis behind its own flag is what makes the forks
independently runnable.

## A note on LD pruning

`--pca_ld_window`, `--pca_ld_step` and `--pca_ld_r2` are defined in
`nextflow.config` but no pruning step is wired in yet. Whether to prune
depends on the analysis: genome-wide structure PCA usually wants LD-pruned
input, while **windowed and locus-specific PCA generally do not** — local LD
is often the signal you are looking for, particularly around inversions. The
parameters are reserved for when a genome-wide PCA module is added.
