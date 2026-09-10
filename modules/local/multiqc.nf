process MULTIQC {
    tag        "cohort"
    label      'process_medium'
    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode

    conda      "bioconda::multiqc=1.25"
    container  "biocontainers/multiqc:1.25--pyhdfd78af_0"

    input:
    path  '*'

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data"       , emit: data
    path "versions.yml"       , emit: versions

    script:
    """
    multiqc --force --title "Cr_align cohort QC" .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """
}
