/*
 * Resolve and, where necessary, build the reference index set.
 *
 * The original pipeline collected `${genome}.*` and destructured the result
 * positionally:
 *     map { bed, fai, pac, sa, amb, genome, ann, bwt, fa, dict -> ... }
 * That binds each index to whatever the glob happened to sort into that slot,
 * so adding or removing a single sidecar file silently mis-assigns every
 * reference input. Here each artefact is resolved by name, and anything
 * missing is built and published for reuse.
 */

include { SAMTOOLS_FAIDX                 } from '../../modules/local/samtools_faidx'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/local/gatk4_createsequencedictionary'
include { BWA_INDEX                      } from '../../modules/local/bwa_index'
include { BUILD_INTERVALS                } from '../../modules/local/build_intervals'

workflow PREPARE_GENOME {

    take:
    fasta          // path
    fai_in         // path or null
    dict_in        // path or null
    bwa_in         // path or null
    chr_regex      // val
    min_length     // val

    main:
    ch_versions = Channel.empty()
    ch_fasta    = Channel.value(file(fasta, checkIfExists: true))

    // ---- .fai ----
    if (fai_in) {
        ch_fai = Channel.value(file(fai_in, checkIfExists: true))
    } else {
        SAMTOOLS_FAIDX(ch_fasta)
        ch_fai      = SAMTOOLS_FAIDX.out.fai.first()
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    // ---- .dict ----
    if (dict_in) {
        ch_dict = Channel.value(file(dict_in, checkIfExists: true))
    } else {
        GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
        ch_dict     = GATK4_CREATESEQUENCEDICTIONARY.out.dict.first()
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }

    // ---- bwa index ----
    if (bwa_in) {
        ch_bwa = Channel.value(file(bwa_in, checkIfExists: true))
    } else {
        BWA_INDEX(ch_fasta)
        ch_bwa      = BWA_INDEX.out.index.first()
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    }

    // ---- scatter intervals ----
    BUILD_INTERVALS(ch_fai, min_length, chr_regex)
    ch_intervals = BUILD_INTERVALS.out.intervals.flatten()

    // Bundle the trio GATK always needs together, so no caller has to
    // reassemble it (and get the order wrong).
    ch_ref = ch_fasta.combine(ch_fai).combine(ch_dict).map { f, i, d -> [f, i, d] }.first()

    emit:
    fasta     = ch_fasta
    fai       = ch_fai
    dict      = ch_dict
    bwa       = ch_bwa
    reference = ch_ref          // [fasta, fai, dict]
    intervals = ch_intervals    // one .interval_list per sequence
    versions  = ch_versions
}
