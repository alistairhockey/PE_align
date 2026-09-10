/*
 * Extract one VCF per named region for locus-specific PCA.
 *
 * The BED's 4th column names each locus; without it, output files would be
 * named by coordinates and the downstream plots would be unlabelled.
 */
process VCF_SUBSET_LOCI {
    tag        "${prefix}"
    label      'process_low'
    publishDir "${params.outdir}/pca/loci", mode: params.publish_dir_mode

    conda      "bioconda::bcftools=1.21"
    container  "quay.io/biocontainers/bcftools:1.21--h8b25389_0"

    input:
    tuple path(vcf), path(tbi)
    path loci_bed
    val  prefix

    output:
    path "*.vcf.gz"    , emit: vcfs
    path "*.vcf.gz.tbi", emit: tbis
    path "loci_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    script:
    """
    printf 'locus\\tregion\\tsites\\tsamples\\n' > loci_summary.tsv

    awk -F'\\t' '!/^#/ && NF>=3 {
        name = (NF>=4 && \$4 != "") ? \$4 : (\$1 "_" \$2 "_" \$3)
        gsub(/[^A-Za-z0-9._-]/, "_", name)
        # BED is 0-based half-open; bcftools -r is 1-based inclusive.
        print name "\\t" \$1 ":" (\$2+1) "-" \$3
    }' ${loci_bed} > regions.tsv

    if [ ! -s regions.tsv ]; then
        echo "ERROR: no usable regions parsed from ${loci_bed}" >&2
        exit 1
    fi

    while IFS=\$'\\t' read -r name region; do
        out="${prefix}.\${name}.vcf.gz"
        bcftools view --threads ${task.cpus} -r "\$region" -Oz -o "\$out" ${vcf}
        bcftools index --tbi "\$out"
        printf '%s\\t%s\\t%s\\t%s\\n' "\$name" "\$region" \\
            "\$(bcftools index -n "\$out")" \\
            "\$(bcftools query -l "\$out" | wc -l)" >> loci_summary.tsv
    done < regions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
