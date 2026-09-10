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

include { NORMALISE_FASTA                } from '../../modules/local/normalise_fasta'
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
    skip_normalise // bool

    main:
    ch_versions = Channel.empty()

    // Normalise line endings and header whitespace before anything indexes
    // the assembly, so every derived artefact agrees on sequence names and
    // coordinates. A clean FASTA passes through untouched.
    //
    // Normalising changes byte offsets, which invalidates any .fai built
    // against the original. So supplying a prebuilt artefact is treated as an
    // assertion that the FASTA is already clean, and normalisation is skipped
    // rather than silently producing a .fai that disagrees with the sequence.
    def has_prebuilt = fai_in || dict_in || bwa_in
    if (skip_normalise || has_prebuilt) {
        if (has_prebuilt && !skip_normalise)
            log.info "PREPARE_GENOME: prebuilt reference artefacts supplied; " +
                     "skipping FASTA normalisation. Ensure the FASTA has LF line " +
                     "endings, or drop --fasta_fai/--fasta_dict/--bwa_index to have " +
                     "the pipeline normalise and rebuild them."
        ch_fasta = Channel.value(file(fasta, checkIfExists: true))
    } else {
        NORMALISE_FASTA(Channel.value(file(fasta, checkIfExists: true)))
        ch_fasta = NORMALISE_FASTA.out.fasta
    }

    // ---- .fai ----
    if (fai_in) {
        ch_fai = Channel.value(file(fai_in, checkIfExists: true))
    } else {
        SAMTOOLS_FAIDX(ch_fasta)
        ch_fai      = SAMTOOLS_FAIDX.out.fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    // ---- .dict ----
    if (dict_in) {
        ch_dict = Channel.value(file(dict_in, checkIfExists: true))
    } else {
        GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
        ch_dict     = GATK4_CREATESEQUENCEDICTIONARY.out.dict
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }

    // ---- bwa index ----
    if (bwa_in) {
        ch_bwa = Channel.value(file(bwa_in, checkIfExists: true))
    } else {
        BWA_INDEX(ch_fasta)
        ch_bwa      = BWA_INDEX.out.index
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    }

    // ---- scatter intervals ----
    // A null --chr_regex means 'every sequence'; '.' matches all names.
    BUILD_INTERVALS(ch_fai, min_length ?: 0, chr_regex ?: '.')
    ch_intervals = BUILD_INTERVALS.out.intervals.flatten()

    // Bundle the trio GATK always needs together, so no caller has to
    // reassemble it (and get the order wrong).
    ch_ref = ch_fasta.combine(ch_fai).combine(ch_dict).map { f, i, d -> [f, i, d] }

    emit:
    fasta     = ch_fasta
    fai       = ch_fai
    dict      = ch_dict
    bwa       = ch_bwa
    reference = ch_ref          // [fasta, fai, dict]
    intervals = ch_intervals    // one .interval_list per sequence
    versions  = ch_versions
}
