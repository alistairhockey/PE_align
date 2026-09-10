process GATK4_GENOTYPEGVCFS {
    tag        "${interval_name}"
    label      'process_high'
    label      'process_high_memory'
    label      'process_long'
    publishDir "${params.outdir}/variants/per_interval", mode: params.publish_dir_mode

    conda      "bioconda::gatk4=4.6.1.0"
    container  "broadinstitute/gatk:4.6.1.0"

    input:
    tuple val(interval_name), path(genomicsdb)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(interval_name), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def avail = task.memory ? (task.memory.giga * 0.7).intValue() : 32
    """
    gatk --java-options "-Xmx${avail}g -XX:-UsePerfData" GenotypeGVCFs \\
        --reference ${fasta} \\
        --variant gendb://${genomicsdb} \\
        --output ${interval_name}.vcf.gz \\
        --create-output-variant-index true \\
        --tmp-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """
}
