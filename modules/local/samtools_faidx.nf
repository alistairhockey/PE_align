process SAMTOOLS_FAIDX {
    tag        "${fasta}"
    label      'process_single'
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode

    conda      "bioconda::samtools=1.21"
    container  "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"

    input:
    path fasta

    output:
    path "*.fai"       , emit: fai
    path "*.gzi"       , emit: gzi, optional: true
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    samtools faidx ${fasta}
    [ -s ${fasta}.fai ] || { echo "ERROR: ${fasta}.fai is missing or empty" >&2; exit 1; }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
