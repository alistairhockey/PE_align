/*
 * Turn an NCBI SRA run table into a resolved run manifest.
 *
 * Runs once per pipeline invocation and is the only step that talks to ENA,
 * so a cohort of any size costs a single request when the run table names its
 * BioProject.
 */
process SRA_RESOLVE {
    tag        "${runtable.name}"
    label      'process_single'
    publishDir "${params.outdir}/metadata", mode: params.publish_dir_mode

    conda      "conda-forge::python=3.11"
    container  "python:3.11-slim"

    input:
    path runtable
    val  layout
    val  assay
    val  organism

    output:
    path "runs.tsv"        , emit: manifest
    path "resolve_summary.json", emit: summary

    script:
    def assay_arg    = assay    ? "--assay '${assay}'"       : ''
    def organism_arg = organism ? "--organism '${organism}'" : ''
    """
    parse_runtable.py ${runtable} \\
        -o runs.tsv \\
        --layout '${layout}' \\
        ${assay_arg} \\
        ${organism_arg} \\
        --summary resolve_summary.json
    """
}
