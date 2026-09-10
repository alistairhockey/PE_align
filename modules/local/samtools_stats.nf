/*
 * Alignment QC. Also feeds the benchmarking report: realised coverage per
 * sample is what turns a wall-clock measurement into a per-Gbp cost model.
 */
process SAMTOOLS_STATS {
    tag        "${meta.id}"
    label      'process_low'
    publishDir "${params.outdir}/qc/samtools", mode: params.publish_dir_mode

    conda      "bioconda::samtools=1.21"
    container  "biocontainers/samtools:1.21--h50ea8bc_0"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("*.stats")   , emit: stats
    tuple val(meta), path("*.flagstat"), emit: flagstat
    tuple val(meta), path("*.coverage"), emit: coverage
    path "versions.yml"                , emit: versions

    script:
    """
    samtools stats    --threads ${task.cpus} --reference ${fasta} ${bam} > ${meta.id}.stats
    samtools flagstat --threads ${task.cpus} ${bam} > ${meta.id}.flagstat
    samtools coverage ${bam} > ${meta.id}.coverage

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
