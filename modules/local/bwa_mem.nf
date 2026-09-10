/*
 * Align one run and emit a coordinate-sorted BAM.
 *
 * The read group is written here and is the single source of truth for sample
 * identity downstream. SM is the *biological sample*, ID is the *run*, so the
 * 77 repeat runs merge into one sample without any later VCF renaming --
 * this is what removes the original pipeline's PullSampleName/RenameVCFs pair.
 */
process BWA_MEM {
    tag        "${meta.id}"
    label      'process_high'

    conda      "bioconda::bwa=0.7.18 bioconda::samtools=1.21"
    container  "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:66ed1b38d280722529bb8a0167b0cf02f8a0b488-0"

    input:
    tuple val(meta), path(reads)
    path  index
    val   fasta_name

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    // Reserve a thread for the sort; bwa and samtools sort run concurrently.
    def align_cpus = Math.max(1, task.cpus - 1)
    def sort_mem   = task.memory ? Math.max(1, ((task.memory.giga * 0.5) / Math.max(1, align_cpus)).intValue()) : 1
    def rg = "@RG\\tID:${meta.run}\\tSM:${meta.sample}\\tPL:ILLUMINA\\tLB:${meta.sample}\\tPU:${meta.run}"
    """
    bwa mem \\
        -t ${align_cpus} \\
        -R '${rg}' \\
        ${index}/${fasta_name} \\
        ${reads} \\
    | samtools sort -@ 1 -m ${sort_mem}G -o ${meta.id}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
