/*
 * Merge one sample's per-interval gVCFs into a single genome-wide gVCF.
 *
 * The scatter produces one gVCF per sample per interval, which is the right
 * shape for GenomicsDBImport but the wrong shape for anything else. Post-hoc
 * population-genetic work -- dxy, windowed divergence, divergence dating --
 * wants one gVCF per sample, retaining the reference-confidence blocks that a
 * joint-called VCF has already collapsed away.
 *
 * MergeVcfs is used rather than bcftools concat because it validates against
 * the sequence dictionary, so a shard from the wrong reference cannot be
 * silently concatenated in.
 */
process GATK4_MERGEVCFS {
    tag        "${sample}"
    label      'process_medium'
    publishDir "${params.outdir}/gvcf", mode: params.publish_dir_mode,
               enabled: params.save_gvcfs

    conda      "bioconda::gatk4=4.6.1.0"
    container  "quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"

    input:
    tuple val(sample), path(gvcfs), path(tbis)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(sample), path("${sample}.g.vcf.gz"), path("${sample}.g.vcf.gz.tbi"), emit: gvcf
    path "versions.yml", emit: versions

    script:
    def avail = task.memory ? (task.memory.giga * 0.7).intValue() : 8
    """
    gatk --java-options "-Xmx${avail}g -XX:-UsePerfData" MergeVcfs \\
        ${gvcfs.collect { "--INPUT ${it}" }.join(' \\\\\n        ')} \\
        --SEQUENCE_DICTIONARY ${dict} \\
        --OUTPUT ${sample}.g.vcf.gz \\
        --CREATE_INDEX true \\
        --TMP_DIR .

    # GATK writes <prefix>.g.vcf.gz.tbi for a bgzipped VCF, but emits
    # <prefix>.g.vcf.tbi in some versions -- normalise.
    [ -f ${sample}.g.vcf.gz.tbi ] || mv ${sample}.g.vcf.tbi ${sample}.g.vcf.gz.tbi

    # A merged gVCF with no records means the shards were empty.
    n=\$(gzip -cd ${sample}.g.vcf.gz | grep -vc '^#' || true)
    if [ "\${n:-0}" -eq 0 ]; then
        echo "ERROR: ${sample}.g.vcf.gz contains no records" >&2
        exit 1
    fi
    echo "${sample}: \$n records across ${gvcfs instanceof List ? gvcfs.size() : 1} intervals"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """
}
