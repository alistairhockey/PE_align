process AlignReads {

    cpus '32'
    memory '64 GB'
    time '24h'


    input:
    tuple val(meta), path(fq)
    tuple path(genome), path(genome_idx)

    output:
    tuple val(meta), path("${meta.srr}.${meta.name}.withDups.bam")

    """
    bwa mem -t ${task.cpus} -R '@RG\\tID:${meta.srr}\\tSM:${meta.name}' ${genome} ${fq} \
    | samtools view -b /dev/stdin \
    | samtools sort /dev/stdin > ${meta.srr}.${meta.name}.withDups.bam
    """
}

process MarkDuplicates {

  memory '64GB'

    publishDir '/scratch/sae001/ahockey/all_ce_align/alignments',  mode: 'copy'

    input:
    tuple val(meta), path(sorted_bam)

    output:
    tuple val(meta), path("${meta.srr}.${meta.name}.bam")

    """
    picard \\
    -Xmx4g \\
    MarkDuplicates \\
   --INPUT ${sorted_bam} \\
   --OUTPUT ${meta.srr}.${meta.name}.bam \\
   --METRICS_FILE metrics.txt \\
   --SORTING_COLLECTION_SIZE_RATIO 0.25 \\
   --READ_NAME_REGEX null

    """
}

process Indexation {

    publishDir '/scratch/sae001/ahockey/all_ce_align/alignments',  mode: 'copy'

    input:
    tuple val(meta), path(markdup_bam)

    output:
    tuple val(meta), path(markdup_bam), path("${markdup_bam}.bai")

    """
    samtools index ${markdup_bam}
    """


}

process MergeBAMs {
  
  memory '64 GB'
  time '24h'

  input:
  tuple val(accession), path(bams), path(bais)

  output:
  tuple val(accession), path("${accession}.bam")

  """
  samtools merge -o ${accession}.merged.bam ${bams} 
  
  samtools sort ${accession}.merged.bam > ${accession}.bam

  """
}

process IdxMerge {

  input:
  tuple val(accession), path(bam)

  output:
  tuple val(accession), path(bam), path("${bam}.bai")

  """
  samtools index ${bam}
  """
}

process HaplotypeCaller {

  memory '64 GB'
  time '24h'

  input:
  tuple val(accession), path(markdupbams), path(markdupbais), val(intName)
  tuple path(genome), path(idx)

  output:
  tuple val(accession), path("${accession}.g.vcf.gz"), val(intName)

  """
  gatk HaplotypeCaller \\
  --reference ${genome} \\
  --emit-ref-confidence GVCF \\
  -I ${markdupbams} \\
  -O ${accession}.g.vcf.gz \\
  -L ${intName} \\  
  """
}

//DO before groupTuple()
process IndexIntGVCFs {

  input:
  tuple val(accession), path(gvcf), val(intName)

  output:
  tuple val(accession), path(gvcf), path("*.tbi"), val(intName)

  """
  gatk IndexFeatureFile \
  -I ${gvcf} 
  """
}



process CombineGVCFs {

  memory '64 GB'
  time '24h'

  input:
  tuple val(accession), path(gvcfs, stageAs: 'input.*.g.vcf.gz'), path(vcfindex, stageAs: 'input.*.g.vcf.gz.tbi'), val(intNames)
  tuple path(genome), path(idx)

  output:
  tuple val(accession), path("*.g.vcf.gz")

  """
  gatk CombineGVCFs \
  -R ${genome} \
  ${gvcfs.collect { "-V $it " }.join()} \
  -O ${accession}.g.vcf.gz 
  """
}

process PullSampleName {

  input:
  tuple val(accession), path(gvcf)

  output:
  tuple val(accession), path(gvcf), path("${accession}.txt")

  """
  bcftools query -l ${gvcf} > ${accession}.txt
  """
  
}

process RenameVCFs{

  input:
  tuple val(accession), path(gvcf), val(sampleName)

  output:
  tuple val(accession), path("${gvcf.simpleName}.sorted.g.vcf.gz")


  """
  picard RenameSampleInVcf \
      INPUT=${gvcf} \
      OUTPUT=${gvcf.simpleName}.sorted.g.vcf.gz \
      OLD_SAMPLE_NAME= ${sampleName} \
      NEW_SAMPLE_NAME=${gvcf.simpleName}
  """
}

process IndexGVCFs {

  publishDir '/scratch/sae001/ahockey/all_ce_align/all_vcfs'

  input:
  tuple val(accession), path(gvcf)

  output:
  tuple val(accession), path(gvcf), path("*.tbi")

  """
  gatk IndexFeatureFile \
  -I ${gvcf} 
  """
}


process ConsolidateGVCFs {
  
  memory '64GB'
  time '72h'

  input:
  tuple val(accessions), path(gvcfs), path(index), val(chrName)

  output:
  tuple val(chrName), path("${chrName}")

  """
  gatk GenomicsDBImport \
  ${gvcfs.collect { "-V $it " }.join()} \
  --intervals ${chrName} \
  --batch-size 50 \
  --genomicsdb-workspace-path ${chrName} 
  """
}
  
process JointCallCohort {

  memory '256GB'
  time '72h'

  publishDir '/scratch/sae001/ahockey/all_ce_align/vcf'

  input:
  tuple val(chrName), path(gvcf_database)
  tuple path(genome), path(idx)


  output:
  tuple val(chrName), path("${chrName}.jointcalledcohort.vcf.gz")

  """
  gatk GenotypeGVCFs \
  --R ${genome} \
  --V gendb://${gvcf_database} \
  --O ${chrName}.jointcalledcohort.vcf.gz
  """
}
