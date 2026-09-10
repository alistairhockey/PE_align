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
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
    END_VERSIONS
    """
}
