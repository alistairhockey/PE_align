/*
 * Gather the per-interval cohort VCFs into one genome-wide call set.
 * Intervals are sorted by their order in the .fai so the output is coordinate
 * sorted without a separate sort pass.
 */
process BCFTOOLS_CONCAT {
    tag        "cohort"
    label      'process_medium'
    publishDir "${params.outdir}/variants", mode: params.publish_dir_mode

    conda      "bioconda::bcftools=1.21"
    container  "quay.io/biocontainers/bcftools:1.21--h8b25389_0"

    input:
    path vcfs
    path tbis
    path fai
    val  prefix

    output:
    tuple path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    # Order the per-interval VCFs to match reference sequence order.
    for v in ${vcfs}; do
        chrom=\$(bcftools view -h "\$v" | grep -v '^##' -m1 >/dev/null; bcftools query -f '%CHROM\\n' "\$v" | head -1)
        idx=\$(grep -n -P "^\${chrom}\\t" ${fai} | cut -d: -f1)
        printf '%s\\t%s\\n' "\${idx:-999999}" "\$v" >> unsorted.tsv
    done
    sort -k1,1n unsorted.tsv | cut -f2 > vcf_order.txt

    bcftools concat \\
        --allow-overlaps \\
        --file-list vcf_order.txt \\
        --threads ${task.cpus} \\
        -Oz -o ${prefix}.vcf.gz
    bcftools index --tbi --threads ${task.cpus} ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
