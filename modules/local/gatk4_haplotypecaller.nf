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

    // Off by default: 1288 shards for the full cohort is ~106 GB, and the
    // merged per-sample gVCFs carry the same information. Enable when you need
    // the scatter units themselves.
    publishDir "${params.outdir}/gvcf/per_interval", mode: params.publish_dir_mode,
               enabled: params.save_interval_gvcfs

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

    # GATK exits 0 after writing a gVCF with no records if the interval file
    # resolved to nothing -- which is exactly what a headerless Picard
    # .interval_list does. An empty gVCF is not a valid result for a whole
    # chromosome, so fail here rather than letting it reach GenomicsDBImport.
    n=\$(gzip -cd ${prefix}.g.vcf.gz | grep -vc '^#' || true)
    if [ "\${n:-0}" -eq 0 ]; then
        echo "ERROR: ${prefix}.g.vcf.gz contains no records." >&2
        echo "The interval file probably resolved to zero bases -- check" >&2
        echo "'Processing N bp from intervals' in the GATK log above." >&2
        exit 1
    fi
    echo "gVCF records: \$n"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """
}
