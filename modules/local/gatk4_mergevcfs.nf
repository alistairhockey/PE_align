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
 *
 * NOTE ON BUILDING ARGUMENT LISTS
 * Do not write  \${files.collect { "--INPUT " + it }.join(" ")}  in a script
 * block. Groovy ends the interpolation at the closure's closing brace, so the
 * remainder of the expression is emitted as literal text and reaches the tool
 * as a stray positional argument. Build such strings in the script: section as
 * a plain variable, or -- as here -- write a list file.
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
    set -euo pipefail

    # Inputs are passed as a '.list' file, one path per line, rather than as
    # repeated --INPUT arguments. A Groovy closure cannot be interpolated into
    # a script block: the closure's closing brace terminates the interpolation
    # early, and the tail of the expression leaks into the command as literal
    # text, which GATK then rejects as an unexpected positional argument.
    # A list file also sidesteps command-line length limits.
    printf '%s\\n' ${gvcfs} > inputs.list
    sort -o inputs.list inputs.list
    echo "merging \$(wc -l < inputs.list) interval gVCFs for ${sample}"

    gatk --java-options "-Xmx${avail}g -XX:-UsePerfData" MergeVcfs \\
        --INPUT inputs.list \\
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
