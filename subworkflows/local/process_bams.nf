/*
 * Turn per-run alignments into analysis-ready per-sample BAMs.
 *
 * Shared by two entry points:
 *   ALIGN_READS   after BWA_MEM, in a normal end-to-end run
 *   --bam_input   when alignments already exist and only the stages from here
 *                 on need to run
 *
 * The second case is not hypothetical: alignment is by far the most expensive
 * stage (~625 GB and ~700 core-hours for 238 runs), so losing it to an
 * unrelated failure downstream -- a full filesystem, a wiped work directory --
 * should not mean realigning.
 */

include { SAMTOOLS_MERGE       } from '../../modules/local/samtools_merge'
include { GATK4_MARKDUPLICATES } from '../../modules/local/gatk4_markduplicates'
include { SAMTOOLS_INDEX       } from '../../modules/local/samtools_index'
include { SAMTOOLS_STATS       } from '../../modules/local/samtools_stats'

workflow PROCESS_BAMS {

    take:
    run_bams     // [meta, bam] one per sequencing run; meta carries .sample
    reference    // [fasta, fai, dict]
    skip_markdup // bool
    skip_qc      // bool

    main:
    ch_versions = Channel.empty()
    ch_qc       = Channel.empty()

    // Group runs into biological samples. 77 of this cohort's 238 runs are
    // second or third runs of a sample, so merging is required, not optional.
    ch_by_sample = run_bams
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

    if (!skip_qc) {
        SAMTOOLS_STATS(ch_bam, reference)
        ch_qc = ch_qc
            .mix(SAMTOOLS_STATS.out.stats.map    { meta, s -> s })
            .mix(SAMTOOLS_STATS.out.flagstat.map { meta, s -> s })
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    }

    emit:
    bam      = ch_bam      // [meta, bam, bai] one per biological sample
    qc       = ch_qc
    versions = ch_versions
}
