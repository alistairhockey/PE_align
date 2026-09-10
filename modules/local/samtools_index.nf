process SAMTOOLS_INDEX {
    tag        "${meta.id}"
    label      'process_single'

    conda      "bioconda::samtools=1.21"
    container  "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("*.bai"), emit: bam
    path "versions.yml"                      , emit: versions

    script:
    """
    samtools index -@ ${task.cpus} ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
