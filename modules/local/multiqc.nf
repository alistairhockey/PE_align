process MULTIQC {
    tag        "cohort"
    label      'process_medium'

    // A QC report is not worth failing a multi-day run over. MultiQC had
    // already written its report successfully when the task was failed on an
    // output-name mismatch, and because the default errorStrategy is 'finish'
    // that stopped every other task from being retried -- including a
    // MarkDuplicates OOM that would otherwise have escalated and succeeded.
    errorStrategy 'ignore'
    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode

    conda      "bioconda::multiqc=1.25"
    container  "quay.io/biocontainers/multiqc:1.25--pyhdfd78af_0"

    input:
    path  '*'

    // MultiQC derives its filenames from --title, so a fixed name here would
    // break whenever the title changes -- which is exactly what happened when
    // the pipeline was renamed and the report became
    // "PE_align-cohort-QC_multiqc_report.html". Globs decouple the two.
    output:
    path "*.html"      , emit: report
    path "*_data"      , emit: data
    path "versions.yml", emit: versions

    script:
    """
    multiqc --force --title "PE_align cohort QC" .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """
}
