/*
 * Produce an analysis-ready SNP matrix for PCA.
 *
 * PCA is sensitive to missingness and to rare variants, and this cohort is
 * median ~8x, so genotype missingness is the dominant artefact. Restricting to
 * biallelic SNPs that pass filters, then applying MAF and missingness
 * thresholds, is what stops PC1 from simply tracking per-sample depth.
 */
process VCF_PREP_PCA {
    tag        "${prefix}"
    label      'process_medium'
    publishDir "${params.outdir}/pca/input", mode: params.publish_dir_mode

    conda      "bioconda::bcftools=1.21"
    container  "quay.io/biocontainers/bcftools:1.21--h8b25389_0"

    input:
    tuple path(vcf), path(tbi)
    val  maf
    val  max_missing
    val  prefix

    output:
    tuple path("${prefix}.pca_input.vcf.gz"), path("${prefix}.pca_input.vcf.gz.tbi"), emit: vcf
    path "${prefix}.pca_input.stats.txt", emit: stats
    path "versions.yml"                 , emit: versions

    script:
    // bcftools expresses the missingness ceiling as an F_MISSING maximum.
    """
    bcftools view \\
        --threads ${task.cpus} \\
        -m2 -M2 -v snps \\
        -f PASS \\
        ${vcf} \\
    | bcftools filter \\
        --threads ${task.cpus} \\
        -e 'F_MISSING > ${max_missing} || MAF < ${maf}' \\
        -Oz -o ${prefix}.pca_input.vcf.gz

    bcftools index --tbi --threads ${task.cpus} ${prefix}.pca_input.vcf.gz

    {
      echo "# PCA input filtering"
      echo "source_vcf:    ${vcf}"
      echo "maf_min:       ${maf}"
      echo "max_missing:   ${max_missing}"
      echo "sites_in:      \$(bcftools index -n ${vcf})"
      echo "sites_out:     \$(bcftools index -n ${prefix}.pca_input.vcf.gz)"
      echo "samples:       \$(bcftools query -l ${prefix}.pca_input.vcf.gz | wc -l)"
    } > ${prefix}.pca_input.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
