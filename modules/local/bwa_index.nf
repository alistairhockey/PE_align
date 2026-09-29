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
    set -euo pipefail
    mkdir bwa
    bwa index -p bwa/${fasta.baseName} ${fasta}

    # `mkdir bwa` succeeds even when indexing does not, so the output glob
    # would match an empty directory and the task would pass. Every downstream
    # BWA_MEM then fails with "fail to locate the index files". Validate that
    # all five index files exist and are non-empty before emitting.
    for ext in amb ann bwt pac sa; do
        f=bwa/${fasta.baseName}.\$ext
        if [ ! -s "\$f" ]; then
            echo "ERROR: bwa index incomplete -- \$f missing or empty" >&2
            ls -la bwa/ >&2
            exit 1
        fi
    done
    echo "bwa index complete:"; ls -la bwa/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p')
    END_VERSIONS
    """
}
