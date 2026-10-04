/*
 * Align every run, merge runs into biological samples, mark duplicates.
 *
 * The merge branch is the reason meta carries `sample` and `run` separately:
 * 75 samples in this cohort have 2 runs and 1 has 3. Single-run samples skip
 * the merge entirely rather than paying for a no-op samtools merge.
 */

include { FASTP                } from '../../modules/local/fastp'
include { FASTQC               } from '../../modules/local/fastqc'
include { BWA_MEM              } from '../../modules/local/bwa_mem'
include { PROCESS_BAMS         } from './process_bams'

workflow ALIGN_READS {

    take:
    reads        // [meta, [fq1, fq2]]
    bwa_index    // path to index dir
    fasta_name   // val: index prefix
    reference    // [fasta, fai, dict]
    trim         // bool
    skip_qc      // bool
    skip_markdup // bool

    main:
    ch_versions = Channel.empty()
    ch_qc       = Channel.empty()

    // ---- Read QC and trimming ----
    if (!skip_qc) {
        FASTQC(reads)
        ch_qc       = ch_qc.mix(FASTQC.out.zip.map { meta, z -> z })
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    if (trim) {
        FASTP(reads)
        ch_to_align = FASTP.out.reads
        ch_qc       = ch_qc.mix(FASTP.out.json.map { meta, j -> j })
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
    } else {
        ch_to_align = reads
    }

    // ---- Align each run ----
    BWA_MEM(ch_to_align, bwa_index, fasta_name)
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    // ---- Per-sample BAMs ----
    // Delegated so the same path can be entered directly from existing
    // alignments via --bam_input. See subworkflows/local/process_bams.nf.
    PROCESS_BAMS(BWA_MEM.out.bam, reference, skip_markdup, skip_qc)
    ch_bam      = PROCESS_BAMS.out.bam
    ch_qc       = ch_qc.mix(PROCESS_BAMS.out.qc)
    ch_versions = ch_versions.mix(PROCESS_BAMS.out.versions)

    emit:
    bam      = ch_bam      // [meta, bam, bai] -- one per biological sample
    qc       = ch_qc
    versions = ch_versions
}
