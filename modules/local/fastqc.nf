process FASTQC {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/qc/fastqc", mode: params.publish_dir_mode

    conda      "bioconda::fastqc=0.12.1"
    container  "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.zip") , emit: zip
    tuple val(meta), path("*.html"), emit: html
    path "versions.yml"            , emit: versions

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """
}
