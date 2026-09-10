/*
 * Merge the runs belonging to one biological sample.
 *
 * 75 of the 161 samples have 2 runs and 1 has 3, so this is a required step,
 * not an optimisation. Single-run samples bypass it in the workflow.
 */
process SAMTOOLS_MERGE {
    tag        "${meta.id}"
    label      'process_medium'

    conda      "bioconda::samtools=1.21"
    container  "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${meta.id}.merged.bam"), emit: bam
    path "versions.yml"                           , emit: versions

    script:
    """
    samtools merge -@ ${task.cpus} -o ${meta.id}.merged.bam ${bams}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
