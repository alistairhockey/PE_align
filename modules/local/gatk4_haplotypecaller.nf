/*
 * Per-sample, per-interval gVCF calling.
 *
 * Scattering here is what makes 161 samples tractable: each task is bounded by
 * one reference sequence rather than the whole genome, and failures retry
 * cheaply. The original passed `-L ${intName}` from splitText(), so interval
 * names arrived with trailing newlines and the command ended on a dangling
 * line continuation.
 */
process GATK4_HAPLOTYPECALLER {
    tag        "${meta.id}:${interval.baseName}"
    label      'process_medium'
    label      'process_long'

    conda      "bioconda::gatk4=4.6.1.0"
    container  "quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"

    input:
    tuple val(meta), path(bam), path(bai), path(interval)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), val("${interval.baseName}"), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), emit: gvcf
    path "versions.yml", emit: versions

    script:
    def avail = task.memory ? (task.memory.giga * 0.8).intValue() : 12
    def prefix = "${meta.id}.${interval.baseName}"
    """
    gatk --java-options "-Xmx${avail}g -XX:-UsePerfData" HaplotypeCaller \\
        --reference ${fasta} \\
        --input ${bam} \\
        --output ${prefix}.g.vcf.gz \\
        --intervals ${interval} \\
        --emit-ref-confidence GVCF \\
        --native-pair-hmm-threads ${task.cpus} \\
        --create-output-variant-index true \\
        --tmp-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """
}
