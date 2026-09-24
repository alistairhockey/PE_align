/*
 * Adapter and quality trimming.
 *
 * --dont_eval_duplication is not optional here. fastp's duplication-rate
 * estimate builds a hash over every read, so its memory grows with read count
 * rather than staying bounded. On Besev_079 (SRR6242343, 82.7 Gbp, ~285M read
 * pairs) that was fatal, and the growth was unmistakable:
 *
 *     32 GB request -> killed at 17.7 GB RSS after  5m23s
 *     64 GB request -> killed at 35.2 GB RSS after 10m44s
 *
 * Twice the memory bought twice the runtime and twice the RSS: unbounded
 * accumulation, not a working set that a larger request would satisfy. Every
 * other run in the cohort peaked near 5 GB, so this is a property of read
 * count, not of the pipeline.
 *
 * Nothing is lost by disabling it. The figure is a QC statistic, and duplicates
 * are actually identified downstream by GATK4_MARKDUPLICATES using alignment
 * positions -- which is both the number that matters and a better estimate than
 * fastp's sequence-identity heuristic.
 */
process FASTP {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/qc/fastp", mode: params.publish_dir_mode, pattern: "*.{json,html}"

    conda      "bioconda::fastp=0.23.4"
    container  "quay.io/biocontainers/fastp:0.23.4--h5f740d0_0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trim.fastq.gz"), emit: reads
    tuple val(meta), path("*.json")         , emit: json
    tuple val(meta), path("*.html")         , emit: html
    path "versions.yml"                     , emit: versions

    script:
    """
    fastp \\
        --in1 ${reads[0]} --in2 ${reads[1]} \\
        --out1 ${meta.id}_1.trim.fastq.gz \\
        --out2 ${meta.id}_2.trim.fastq.gz \\
        --detect_adapter_for_pe \\
        --dont_eval_duplication \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
    END_VERSIONS
    """
}
