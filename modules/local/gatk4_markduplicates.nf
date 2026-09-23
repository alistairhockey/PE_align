/*
 * Mark optical/PCR duplicates.
 *
 * Java heap is derived from task.memory rather than the original's fixed
 * -Xmx4g inside a 64 GB request, which both wasted the allocation and left
 * Picard's sorting collection to spill to disk unnecessarily.
 */
process GATK4_MARKDUPLICATES {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/alignments",   mode: params.publish_dir_mode, pattern: "*.{bam,bai}"
    publishDir "${params.outdir}/qc/markdup",   mode: params.publish_dir_mode, pattern: "*.metrics.txt"

    conda      "bioconda::gatk4=4.6.1.0"
    container  "quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    tuple val(meta), path("*.metrics.txt")                             , emit: metrics
    path "versions.yml"                                                , emit: versions

    script:
    // 70%, not 80%. Picard's sorting collection, the I/O buffers and the JVM's
    // own native allocation all live outside the heap, and at 64 GB an -Xmx51g
    // heap left too little for them: Bari1_092 (~73x coverage, 340M+ records)
    // was killed by the cgroup mid-write with exit 247.
    def avail = task.memory ? (task.memory.giga * 0.7).intValue() : 8
    """
    gatk --java-options "-Xmx${avail}g -XX:-UsePerfData" MarkDuplicates \\
        --INPUT ${bam} \\
        --OUTPUT ${meta.id}.bam \\
        --METRICS_FILE ${meta.id}.metrics.txt \\
        --CREATE_INDEX true \\
        --SORTING_COLLECTION_SIZE_RATIO 0.25 \\
        --READ_NAME_REGEX null \\
        --VALIDATION_STRINGENCY LENIENT \\
        --TMP_DIR .

    # GATK writes <prefix>.bai; downstream tools and IGV expect <prefix>.bam.bai
    mv ${meta.id}.bai ${meta.id}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """
}
