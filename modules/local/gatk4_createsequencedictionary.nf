process GATK4_CREATESEQUENCEDICTIONARY {
    tag        "${fasta}"
    label      'process_low'
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode

    conda      "bioconda::gatk4=4.6.1.0"
    container  "quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"

    input:
    path fasta

    output:
    path "*.dict"      , emit: dict
    path "versions.yml", emit: versions

    script:
    // Give the JVM most of the request but leave headroom for native/off-heap.
    def avail = task.memory ? (task.memory.giga * 0.8).intValue() : 4
    """
    set -euo pipefail
    gatk --java-options "-Xmx${avail}g" CreateSequenceDictionary \\
        --REFERENCE ${fasta} \\
        --OUTPUT ${fasta.baseName.replaceAll(/\.(fa|fasta|fna)$/, "")}.dict

    d=${fasta.baseName.replaceAll(/\.(fa|fasta|fna)$/, "")}.dict
    grep -q '^@SQ' "\$d" || { echo "ERROR: \$d has no @SQ records" >&2; exit 1; }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """
}
