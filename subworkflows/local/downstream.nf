/*
 * Analysis-ready outputs for downstream population-genetic work.
 *
 * These are deliberately *fork points*, not analyses. The pipeline produces a
 * filtered SNP matrix and the per-chromosome / per-locus slices that windowed
 * and locus-specific PCA consume; the PCA scripts themselves stay external so
 * they can be iterated on without rerunning variant calling.
 *
 *   --run_winpca      per-chromosome VCFs      -> results/pca/winpca/
 *   --run_locus_pca   per-locus VCFs (--loci_bed) -> results/pca/loci/
 */

include { VCF_PREP_PCA     } from '../../modules/local/vcf_prep_pca'
include { VCF_SPLIT_CHROM  } from '../../modules/local/vcf_split_chrom'
include { VCF_SUBSET_LOCI  } from '../../modules/local/vcf_subset_loci'

workflow DOWNSTREAM {

    take:
    vcf            // [vcf, tbi] -- the filtered cohort call set
    fai            // path
    chr_regex      // val
    run_winpca     // bool
    run_locus_pca  // bool
    loci_bed       // path or null
    prefix         // val

    main:
    ch_versions = Channel.empty()

    VCF_PREP_PCA(vcf, params.pca_maf, params.pca_max_missing, prefix)
    ch_versions = ch_versions.mix(VCF_PREP_PCA.out.versions)

    ch_winpca = Channel.empty()
    if (run_winpca) {
        VCF_SPLIT_CHROM(VCF_PREP_PCA.out.vcf, fai, chr_regex ?: '.', prefix)
        ch_winpca   = VCF_SPLIT_CHROM.out.vcfs
        ch_versions = ch_versions.mix(VCF_SPLIT_CHROM.out.versions)
    }

    ch_loci = Channel.empty()
    if (run_locus_pca) {
        if (!loci_bed)
            error "--run_locus_pca requires --loci_bed <regions.bed> (4th column = locus name)"
        VCF_SUBSET_LOCI(VCF_PREP_PCA.out.vcf, file(loci_bed, checkIfExists: true), prefix)
        ch_loci     = VCF_SUBSET_LOCI.out.vcfs
        ch_versions = ch_versions.mix(VCF_SUBSET_LOCI.out.versions)
    }

    emit:
    pca_input = VCF_PREP_PCA.out.vcf
    winpca    = ch_winpca
    loci      = ch_loci
    versions  = ch_versions
}
