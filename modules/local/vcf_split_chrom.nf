/*
 * Split the PCA input into per-chromosome VCFs.
 *
 * winPCA slides a window along one sequence at a time, so per-chromosome
 * inputs let the downstream analysis parallelise and keep memory bounded.
 */
process VCF_SPLIT_CHROM {
    tag        "${prefix}"
    label      'process_medium'
    publishDir "${params.outdir}/pca/winpca", mode: params.publish_dir_mode

    conda      "bioconda::bcftools=1.21"
    container  "biocontainers/bcftools:1.21--h8b25389_0"

    input:
    tuple path(vcf), path(tbi)
    path fai
    val  chr_regex
    val  prefix

    output:
    path "*.vcf.gz"    , emit: vcfs
    path "*.vcf.gz.tbi", emit: tbis
    path "chrom_list.txt", emit: chroms
    path "versions.yml", emit: versions

    script:
    def re = chr_regex ?: '.'
    """
    awk -F'\\t' '\$1 ~ /${re}/ {print \$1}' ${fai} > chrom_list.txt

    while read -r chrom; do
        safe=\$(echo "\$chrom" | tr -c 'A-Za-z0-9._-' '_')
        bcftools view --threads ${task.cpus} -r "\$chrom" -Oz \\
            -o "${prefix}.\${safe}.vcf.gz" ${vcf}
        bcftools index --tbi "${prefix}.\${safe}.vcf.gz"
    done < chrom_list.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
