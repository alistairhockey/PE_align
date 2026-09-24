/*
 * Scatter-gather joint genotyping.
 *
 *   HaplotypeCaller   per sample x interval    -> gVCF
 *   GenomicsDBImport  per interval, all samples-> workspace
 *   GenotypeGVCFs     per interval             -> cohort VCF shard
 *   bcftools concat   all shards               -> genome-wide cohort VCF
 *
 * The original combined gVCFs per *sample* across intervals before importing,
 * which serialised work that the scatter had just parallelised and produced
 * no artefact the cohort join needed. Sample renaming is also gone: the RG SM
 * tag set at alignment already carries the correct name.
 */

include { GATK4_HAPLOTYPECALLER   } from '../../modules/local/gatk4_haplotypecaller'
include { GATK4_MERGEVCFS         } from '../../modules/local/gatk4_mergevcfs'
include { GATK4_GENOMICSDBIMPORT  } from '../../modules/local/gatk4_genomicsdbimport'
include { GATK4_GENOTYPEGVCFS     } from '../../modules/local/gatk4_genotypegvcfs'
include { BCFTOOLS_CONCAT         } from '../../modules/local/bcftools_concat'
include { GATK4_VARIANTFILTRATION } from '../../modules/local/gatk4_variantfiltration'

workflow CALL_VARIANTS {

    take:
    bam            // [meta, bam, bai]
    reference      // [fasta, fai, dict]
    fai            // path
    intervals      // one .interval_list per sequence
    skip_filtering // bool
    prefix         // val

    main:
    ch_versions = Channel.empty()

    // Full cross of samples x intervals.
    ch_hc_input = bam.combine(intervals)

    GATK4_HAPLOTYPECALLER(ch_hc_input, reference)
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

    // Per-sample genome-wide gVCFs, kept for post-hoc population genetics
    // (dxy, windowed divergence, divergence dating). These retain the
    // reference-confidence blocks that joint calling collapses away, so they
    // cannot be reconstructed from the cohort VCF.
    ch_by_sample = GATK4_HAPLOTYPECALLER.out.gvcf
        .map { meta, interval_name, gvcf, tbi -> [ meta.id, gvcf, tbi ] }
        .groupTuple()

    GATK4_MERGEVCFS(ch_by_sample, reference)
    ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions.first())

    // Re-key by interval and gather every sample's shard for that interval.
    ch_interval_files = intervals.map { iv -> [ iv.baseName, iv ] }

    ch_gendb_input = GATK4_HAPLOTYPECALLER.out.gvcf
        .map { meta, interval_name, gvcf, tbi -> [ interval_name, gvcf, tbi ] }
        .groupTuple()
        .join(ch_interval_files)

    GATK4_GENOMICSDBIMPORT(ch_gendb_input, reference)
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions.first())

    GATK4_GENOTYPEGVCFS(GATK4_GENOMICSDBIMPORT.out.genomicsdb, reference)
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions.first())

    BCFTOOLS_CONCAT(
        GATK4_GENOTYPEGVCFS.out.vcf.map { n, v, t -> v }.collect(),
        GATK4_GENOTYPEGVCFS.out.vcf.map { n, v, t -> t }.collect(),
        fai,
        prefix
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    ch_raw = BCFTOOLS_CONCAT.out.vcf

    if (!skip_filtering) {
        GATK4_VARIANTFILTRATION(
            ch_raw, reference,
            params.snp_filter_expression,
            params.indel_filter_expression,
            prefix
        )
        ch_final    = GATK4_VARIANTFILTRATION.out.vcf
        ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions)
    } else {
        ch_final = ch_raw
    }

    emit:
    gvcf     = GATK4_MERGEVCFS.out.gvcf     // [sample, gvcf, tbi] per sample
    raw      = ch_raw
    vcf      = ch_final
    versions = ch_versions
}
