/*
 * Consolidate every sample's gVCF for one interval into a GenomicsDB workspace.
 *
 * This is the cohort-wide join, so it is the single largest memory consumer in
 * the pipeline and the step most likely to dictate the Pawsey request.
 * --batch-size bounds the number of readers held open at once; without it,
 * 161 samples exhaust file handles before they exhaust memory.
 */
process GATK4_GENOMICSDBIMPORT {
    tag        "${interval_name}"
    label      'process_high'
    label      'process_high_memory'
    label      'process_long'

    conda      "bioconda::gatk4=4.6.1.0"
    container  "broadinstitute/gatk:4.6.1.0"

    input:
    tuple val(interval_name), path(gvcfs), path(tbis), path(interval)

    output:
    tuple val(interval_name), path("${interval_name}_gendb"), emit: genomicsdb
    path "versions.yml", emit: versions

    script:
    def avail = task.memory ? (task.memory.giga * 0.7).intValue() : 32
    """
    # A sample map keeps the command line bounded regardless of cohort size.
    for f in ${gvcfs}; do
        sample=\$(gzip -cd "\$f" | grep -m1 '^#CHROM' | cut -f10-)
        printf '%s\\t%s\\n' "\$sample" "\$f" >> sample_map.tsv
    done

    gatk --java-options "-Xmx${avail}g -Xms${Math.max(1, (avail/2).intValue())}g -XX:-UsePerfData" \\
        GenomicsDBImport \\
        --sample-name-map sample_map.tsv \\
        --genomicsdb-workspace-path ${interval_name}_gendb \\
        --intervals ${interval} \\
        --batch-size 50 \\
        --reader-threads ${task.cpus} \\
        --merge-input-intervals \\
        --tmp-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """
}
