process BWA_INDEX {
    tag        "${fasta.baseName}"
    label      'process_high'
    publishDir "${params.outdir}/reference/bwa", mode: params.publish_dir_mode

    conda      "bioconda::bwa=0.7.18"
    container  "quay.io/biocontainers/bwa:0.7.18--he4a0461_1"

    input:
    path fasta

    output:
    path "bwa"         , emit: index
    path "versions.yml", emit: versions

    script:
    """
    mkdir bwa
    bwa index -p bwa/${fasta.baseName} ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p')
    END_VERSIONS
    """
}
