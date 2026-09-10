/*
 * Download one run's paired FASTQs from ENA and verify both checksums.
 *
 * `storeDir` makes this idempotent across runs: if both files are already
 * present in the reads directory the task never executes, so re-running the
 * pipeline never re-downloads. A file that fails MD5 is deleted and the task
 * fails, so a corrupt download is retried by the normal errorStrategy rather
 * than silently poisoning the alignment.
 */
process SRA_FETCH {
    tag      "${meta.run}"
    label    'process_single'
    storeDir "${params.reads_dir ?: "${params.outdir}/reads"}"

    // Downloads are network-bound, not CPU-bound; too many at once makes ENA
    // throttle everything. maxForks is set from params.max_download_jobs.
    maxForks { params.max_download_jobs as int }
    errorStrategy 'retry'
    maxRetries    3

    conda     "conda-forge::python=3.11"
    container "python:3.11-slim"

    input:
    tuple val(meta), val(url_1), val(md5_1), val(url_2), val(md5_2)

    output:
    tuple val(meta), path("${meta.run}_${meta.sample}_1.fastq.gz"),
                     path("${meta.run}_${meta.sample}_2.fastq.gz"), emit: reads

    script:
    def o1 = "${meta.run}_${meta.sample}_1.fastq.gz"
    def o2 = "${meta.run}_${meta.sample}_2.fastq.gz"
    """
    fetch_run.py '${url_1}' '${md5_1}' '${o1}'
    fetch_run.py '${url_2}' '${md5_2}' '${o2}'
    """
}
