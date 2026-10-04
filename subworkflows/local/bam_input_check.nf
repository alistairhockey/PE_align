/*
 * Parse a samplesheet of existing per-run alignments.
 *
 * Columns: sample,run,bam
 *   sample  biological sample -- runs sharing it are merged
 *   run     sequencing run, unique within a sample
 *   bam     per-run BAM, coordinate sorted, as produced by BWA_MEM
 *
 * The BAMs must carry correct @RG SM tags, since everything downstream reads
 * sample identity from the BAM header rather than from this sheet. BAMs written
 * by this pipeline's BWA_MEM always do.
 */

workflow BAM_INPUT_CHECK {

    take:
    samplesheet   // path
    subset        // int or null

    main:
    ch_bams = Channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            [ 'sample', 'run', 'bam' ].each { col ->
                if (!row.containsKey(col))
                    error "--bam_input sheet is missing column '${col}'. Found: ${row.keySet().join(', ')}"
            }
            if (!row.sample?.trim()) error "--bam_input row has an empty 'sample': ${row}"
            if (!row.run?.trim())    error "--bam_input row has an empty 'run': ${row}"
            def meta = [
                id        : "${row.sample}_${row.run}".toString(),
                sample    : row.sample.trim(),
                run       : row.run.trim(),
                single_end: false
            ]
            [ meta, file(row.bam, checkIfExists: true) ]
        }

    if (subset) {
        ch_bams = ch_bams
            .toList()
            .flatMap { rows ->
                def keep = rows.collect { it[0].sample }.unique().sort().take(subset as int)
                rows.findAll { keep.contains(it[0].sample) }
            }
    }

    ch_bams
        .map { meta, bam -> meta.id }
        .toList()
        .map { ids ->
            def dups = ids.countBy { it }.findAll { k, v -> v > 1 }.keySet()
            if (dups) error "Duplicated sample+run in --bam_input: ${dups.join(', ')}"
            ids.size()
        }
        .subscribe { n -> log.info "BAM_INPUT_CHECK: ${n} per-run BAM(s) accepted" }

    emit:
    bams = ch_bams   // [meta, bam]
}
