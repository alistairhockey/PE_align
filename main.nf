#!/usr/bin/env nextflow
/*
========================================================================================
    PE_align
========================================================================================
    Paired-end short-read alignment and joint variant calling, instrumented
    for resource benchmarking.

    Developed for a 238-run Cicer reticulatum WGS cohort against the
    PBA_HatTrick chickpea assembly, but nothing in it is specific to that
    cohort or that reference.

    Two inputs and go:
        --sra_metadata   an NCBI SRA run table (SraRunTable.txt)
        --fasta          a genome assembly

    Everything else -- FASTQ downloads, .fai, .dict, bwa index, scatter
    intervals, the samplesheet -- is derived.

    https://github.com/alistairhockey/PE_align
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { FETCH_READS    } from './subworkflows/local/fetch_reads'
include { INPUT_CHECK    } from './subworkflows/local/input_check'
include { PREPARE_GENOME } from './subworkflows/local/prepare_genome'
include { ALIGN_READS    } from './subworkflows/local/align_reads'
include { CALL_VARIANTS  } from './subworkflows/local/call_variants'
include { DOWNSTREAM     } from './subworkflows/local/downstream'
include { MULTIQC        } from './modules/local/multiqc'

/*
========================================================================================
    HELP
========================================================================================
*/

def helpMessage() {
    log.info """
    ==========================================================================
     PE_align  v${workflow.manifest.version}
    ==========================================================================

    Typical use:

      nextflow run . -profile uwa,apptainer \\
          --sra_metadata SraRunTable.txt \\
          --fasta /path/to/assembly.fna \\
          --reads_dir /group/peg/cicer/cret/reads \\
          --outdir results

    REQUIRED
      --fasta               Genome assembly FASTA. Index files, sequence
                            dictionary and scatter intervals are built from it
                            and published to <outdir>/reference/.
      One of:
      --sra_metadata        NCBI SRA run table; runs are resolved against ENA
                            and downloaded.
      --input               Samplesheet CSV: sample,run,fastq_1,fastq_2

    READS
      --reads_dir           Reuse FASTQs already present here. Files are named
                            <RUN>_<SAMPLE>_{1,2}.fastq.gz
      --download_reads      Download anything missing (default: true)
      --max_download_jobs   Concurrent downloads (default: 8)
      --sra_layout          Keep runs with this layout (default: PAIRED)
      --sra_assay           Keep runs with this assay type, e.g. WGS
      --sra_organism        Keep runs whose ENA organism matches

    REFERENCE
      --fasta_fai --fasta_dict --bwa_index    Reuse prebuilt artefacts
      --chr_regex           Only call on sequences matching this regex
      --intervals_min_length   Skip sequences shorter than this (bp)

    STAGES
      --trim_reads          fastp adapter/quality trimming (default: true)
      --skip_qc --skip_markduplicates --skip_variant_calling
      --skip_filtering --skip_multiqc

    DOWNSTREAM FORKS
      --run_winpca          Per-chromosome VCFs for windowed PCA
      --run_locus_pca       Per-locus VCFs; needs --loci_bed
      --loci_bed            BED of regions; 4th column names each locus
      --pca_maf --pca_max_missing

    BENCHMARKING
      --benchmark_label     Tag for this run's trace/report/timeline files
      --bench_subset        Use only the first N samples
      Profiles bench_2 / bench_8 / bench_24 set --bench_subset for scaling runs.

    PROFILES
      uwa, setonix, slurm, standard          execution environment
      apptainer, singularity, docker, conda  software provisioning
    ==========================================================================
    """.stripIndent()
}

/*
========================================================================================
    PARAMETER VALIDATION
========================================================================================
*/

if (params.help) { helpMessage(); exit 0 }

def errors = []

if (!params.fasta)
    errors << "--fasta is required (genome assembly FASTA)."
else if (!file(params.fasta).exists())
    errors << "--fasta not found: ${params.fasta}"

if (!params.sra_metadata && !params.input)
    errors << "Provide either --sra_metadata <SraRunTable.txt> or --input <samplesheet.csv>."
if (params.sra_metadata && params.input)
    errors << "--sra_metadata and --input are mutually exclusive; pick one."
if (params.sra_metadata && !file(params.sra_metadata).exists())
    errors << "--sra_metadata not found: ${params.sra_metadata}"
if (params.input && !file(params.input).exists())
    errors << "--input not found: ${params.input}"
if (params.run_locus_pca && !params.loci_bed)
    errors << "--run_locus_pca requires --loci_bed."

if (errors) {
    log.error "Parameter validation failed:\n  - " + errors.join("\n  - ") +
              "\n\nRun with --help for usage."
    exit 1
}

/*
========================================================================================
    RUN SUMMARY
========================================================================================
*/

def subset = params.containsKey('bench_subset') ? params.bench_subset : null

log.info """
==========================================================================
 PE_align v${workflow.manifest.version}
==========================================================================
 input        : ${params.sra_metadata ?: params.input}
 reference    : ${params.fasta}
 reads_dir    : ${params.reads_dir ?: "${params.outdir}/reads"}
 outdir       : ${params.outdir}
 benchmarks   : ${params.tracedir}
 profile      : ${workflow.profile}
 subset       : ${subset ?: 'none (full cohort)'}
 downstream   : winPCA=${params.run_winpca} locusPCA=${params.run_locus_pca}
==========================================================================
""".stripIndent()

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    ch_versions = Channel.empty()

    // ---- Reads ----
    if (params.sra_metadata) {
        FETCH_READS(
            file(params.sra_metadata),
            params.reads_dir,
            params.download_reads,
            params.sra_layout   ?: 'PAIRED',
            params.sra_assay,
            params.sra_organism,
            subset
        )
        ch_reads = FETCH_READS.out.reads
    } else {
        INPUT_CHECK(file(params.input), subset)
        ch_reads = INPUT_CHECK.out.reads
    }

    // ---- Reference ----
    PREPARE_GENOME(
        params.fasta,
        params.fasta_fai,
        params.fasta_dict,
        params.bwa_index,
        params.chr_regex,
        params.intervals_min_length,
        params.skip_fasta_normalisation
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // bwa index files are named after the FASTA, minus any .gz
    def fasta_name = file(params.fasta).name.replaceAll(/\.gz$/, '')
                        .replaceAll(/\.(fa|fasta|fna)$/, '')
    def prefix     = params.benchmark_label ?: 'cohort'

    // ---- Align ----
    ALIGN_READS(
        ch_reads,
        PREPARE_GENOME.out.bwa,
        fasta_name,
        PREPARE_GENOME.out.reference,
        params.trim_reads,
        params.skip_qc,
        params.skip_markduplicates
    )
    ch_versions = ch_versions.mix(ALIGN_READS.out.versions)

    // ---- Call ----
    if (!params.skip_variant_calling) {

        CALL_VARIANTS(
            ALIGN_READS.out.bam,
            PREPARE_GENOME.out.reference,
            PREPARE_GENOME.out.fai,
            PREPARE_GENOME.out.intervals,
            params.skip_filtering,
            prefix
        )
        ch_versions = ch_versions.mix(CALL_VARIANTS.out.versions)

        // ---- Downstream forks ----
        if (params.run_winpca || params.run_locus_pca) {
            DOWNSTREAM(
                CALL_VARIANTS.out.vcf,
                PREPARE_GENOME.out.fai,
                params.chr_regex,
                params.run_winpca,
                params.run_locus_pca,
                params.loci_bed,
                prefix
            )
            ch_versions = ch_versions.mix(DOWNSTREAM.out.versions)
        }
    }

    // ---- Reporting ----
    if (!params.skip_multiqc) {
        MULTIQC(ALIGN_READS.out.qc.collect().ifEmpty([]))
    }
}

/*
========================================================================================
    COMPLETION
========================================================================================
*/

workflow.onComplete {
    def status = workflow.success ? 'COMPLETED' : 'FAILED'
    log.info """
==========================================================================
 PE_align ${status}
--------------------------------------------------------------------------
 duration     : ${workflow.duration}
 CPU hours    : ${workflow.stats?.computeTimeFmt ?: 'n/a'}
 tasks ok     : ${workflow.stats?.succeededCount ?: 0}
 tasks failed : ${workflow.stats?.failedCount ?: 0}
 tasks cached : ${workflow.stats?.cachedCount ?: 0}
 results      : ${params.outdir}
 benchmarks   : ${params.tracedir}

 Next: summarise resource use for the allocation case with
   bin/summarise_benchmark.py ${params.tracedir}/trace-*.txt
==========================================================================
""".stripIndent()
}
