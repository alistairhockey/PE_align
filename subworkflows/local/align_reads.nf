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
include { SAMTOOLS_MERGE       } from '../../modules/local/samtools_merge'
include { GATK4_MARKDUPLICATES } from '../../modules/local/gatk4_markduplicates'
include { SAMTOOLS_INDEX       } from '../../modules/local/samtools_index'
include { SAMTOOLS_STATS       } from '../../modules/local/samtools_stats'

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

    // ---- Group runs into samples ----
    // Rewrite meta so identity is the sample, not the sample+run pair.
    ch_by_sample = BWA_MEM.out.bam
        .map { meta, bam -> [ [id: meta.sample, sample: meta.sample], bam ] }
        .groupTuple()
        .branch { meta, bams ->
            single:   bams.size() == 1
                return [ meta, bams[0] ]
            multiple: bams.size() > 1
                return [ meta, bams ]
        }

    SAMTOOLS_MERGE(ch_by_sample.multiple)
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    ch_sample_bam = ch_by_sample.single.mix(SAMTOOLS_MERGE.out.bam)

    // ---- Duplicates ----
    if (!skip_markdup) {
        GATK4_MARKDUPLICATES(ch_sample_bam)
        ch_bam      = GATK4_MARKDUPLICATES.out.bam
        ch_qc       = ch_qc.mix(GATK4_MARKDUPLICATES.out.metrics.map { meta, m -> m })
        ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())
    } else {
        SAMTOOLS_INDEX(ch_sample_bam)
        ch_bam      = SAMTOOLS_INDEX.out.bam
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    }

    // ---- Alignment QC ----
    if (!skip_qc) {
        SAMTOOLS_STATS(ch_bam, reference)
        ch_qc = ch_qc
            .mix(SAMTOOLS_STATS.out.stats.map    { meta, s -> s })
            .mix(SAMTOOLS_STATS.out.flagstat.map { meta, s -> s })
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    }

    emit:
    bam      = ch_bam      // [meta, bam, bai] -- one per biological sample
    qc       = ch_qc
    versions = ch_versions
}
