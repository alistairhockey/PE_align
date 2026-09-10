
include { AlignReads } from params.modules_path
include { MarkDuplicates } from params.modules_path
include { Indexation } from params.modules_path
include { MergeBAMs } from params.modules_path
include { IdxMerge } from params.modules_path
include { HaplotypeCaller } from params.modules_path
include { IndexIntGVCFs } from params.modules_path
include { CombineGVCFs } from params.modules_path
include { PullSampleName } from params.modules_path
include { RenameVCFs } from params.modules_path
include { IndexGVCFs } from params.modules_path
include { ConsolidateGVCFs } from params.modules_path
include { JointCallCohort } from params.modules_path


workflow {

    Channel.fromFilePairs(params.pe_reads, size: params.singleEnd ? 1 : 2 ) 
    // // Channel.fromPath(params.reads)
    | map {it, path ->
        def idParts = it.split("_Wild_Cicer_Collection_")
         [[srr: idParts[0], name: idParts[1]], path]}
    | set { pe_reads }

    Channel.fromPath(params.all_reads) 
       | filter {!(it.name =~ /.*_[12]\.fastq\.gz$/) } 
    | map { it -> 
            def idParts = it.simpleName.split("_Wild_Cicer_Collection_")
            [[srr: idParts[0], name: idParts[1]], it]}     // | flatten
    | set { se_reads }
    
    pe_reads
    | concat(se_reads)
    | set { samples }


    Channel.fromPath("${params.genome}.*") 
    | collect 
    | map {bed, fai, pac, sa, amb, genome, ann, bwt, fa, dict -> [[genome], [fai, pac, sa, amb, ann, bwt, dict]]}
    | first
    | set {genome}

    AlignReads (samples, genome)
    | MarkDuplicates
    | Indexation
    | set {bams}

    bams
    | map { [it[0].name, it[1..2]] }  // Transform the data into tuples
    | groupTuple()  // Group the tuples by name
    | filter { it[1].size() > 1 } 
    | map { id, files ->
        (bams, bais) = files.transpose()
        [id, bams, bais]
        }
    | set {to_be_merged}

    to_be_merged
    | MergeBAMs
    | IdxMerge
    | set {merged_bams}

    bams
    | map { [it[0].name, it[1..2]] }  // Transform the data into tuples
    | groupTuple()  // Group the tuples by name
    | filter { it[1].size() == 1 } 
    | map { id, files -> [id, files[0][0], files[0][1]]}
    | concat(merged_bams)
    | set {all_bams}
    
    

    Channel.fromPath(params.intervals) 
    | splitText()
    | set {chromosomeIntervals}

    all_bams
    | combine (chromosomeIntervals)
    | set {alignments}

    HaplotypeCaller(alignments, genome)
    | IndexIntGVCFs
    | groupTuple()
    | set {gvcf}

    Channel.fromPath("${params.genome}.fna.fai")
    | splitCsv(header:["chr", "stopbase", "a", "b", "c"], sep: "\t") 
    | filter { row -> row.chr =~ /^cicec.S2Drd065.gnm1.chr*/ }
    | map { row -> row.chr }
    | set { chromosomeNames }

    CombineGVCFs(gvcf, genome)
    | PullSampleName
    | map {id, vcf, file ->[id, vcf, file.text.trim()] } 
    | RenameVCFs 
    | IndexGVCFs 
    | combine(chromosomeNames) 
    | map { id, vcf, chrName ->
        def clean_chrName = chrName.replace('|', '_') 
        [ id, vcf, clean_chrName ]
    }
    | groupTuple(by: 3) | view
    // | ConsolidateGVCFs 
    // | set {consolidate}

    // JointCallCohort(consolidate, genome)



}