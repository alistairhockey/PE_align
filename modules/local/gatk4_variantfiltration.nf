/*
 * GATK hard filtering, applied separately to SNPs and indels because the
 * recommended thresholds differ. VQSR is not used: it needs a truth set that
 * does not exist for wild Cicer.
 */
process GATK4_VARIANTFILTRATION {
    tag        "cohort"
    label      'process_medium'
    label      'process_long'
    publishDir "${params.outdir}/variants", mode: params.publish_dir_mode

    conda      "bioconda::gatk4=4.6.1.0"
    container  "broadinstitute/gatk:4.6.1.0"

    input:
    tuple path(vcf), path(tbi)
    tuple path(fasta), path(fai), path(dict)
    val  snp_expr
    val  indel_expr
    val  prefix

    output:
    tuple path("${prefix}.filtered.vcf.gz"), path("${prefix}.filtered.vcf.gz.tbi"), emit: vcf
    path "*.filter_summary.txt", emit: summary
    path "versions.yml"        , emit: versions

    script:
    def avail = task.memory ? (task.memory.giga * 0.8).intValue() : 16
    """
    gatk --java-options "-Xmx${avail}g" SelectVariants \\
        -R ${fasta} -V ${vcf} --select-type-to-include SNP -O snps.vcf.gz

    gatk --java-options "-Xmx${avail}g" SelectVariants \\
        -R ${fasta} -V ${vcf} --select-type-to-include INDEL \\
        --select-type-to-include MIXED -O indels.vcf.gz

    gatk --java-options "-Xmx${avail}g" VariantFiltration \\
        -R ${fasta} -V snps.vcf.gz \\
        --filter-name "SNP_HARD_FILTER" --filter-expression "${snp_expr}" \\
        -O snps.filt.vcf.gz

    gatk --java-options "-Xmx${avail}g" VariantFiltration \\
        -R ${fasta} -V indels.vcf.gz \\
        --filter-name "INDEL_HARD_FILTER" --filter-expression "${indel_expr}" \\
        -O indels.filt.vcf.gz

    gatk --java-options "-Xmx${avail}g" MergeVcfs \\
        -I snps.filt.vcf.gz -I indels.filt.vcf.gz \\
        -O ${prefix}.filtered.vcf.gz

    {
      echo "# Variant filtering summary"
      echo "snp_expression:   ${snp_expr}"
      echo "indel_expression: ${indel_expr}"
      echo
      printf 'total\\t%s\\n'  "\$(gatk CountVariants -V ${prefix}.filtered.vcf.gz 2>/dev/null | tail -1)"
    } > ${prefix}.filter_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """
}
