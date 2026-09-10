/*
 * Parse and validate the samplesheet.
 *
 * Samplesheet columns: sample,run,fastq_1,fastq_2
 *   sample  biological sample -- becomes the RG SM tag and the BAM/VCF name
 *   run     sequencing run    -- becomes the RG ID; unique within a sample
 *
 * A sample with several runs appears as several rows. This is the norm for
 * this cohort: 238 runs across 161 samples.
 */

workflow INPUT_CHECK {

    take:
    samplesheet   // path
    subset        // int or null: keep only the first N samples (benchmarking)

    main:
    ch_reads = Channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            [ 'sample', 'run', 'fastq_1', 'fastq_2' ].each { col ->
                if (!row.containsKey(col))
                    error "Samplesheet is missing required column '${col}'. Found: ${row.keySet().join(', ')}"
            }
            if (!row.sample?.trim()) error "Samplesheet has a row with an empty 'sample' value: ${row}"
            if (!row.run?.trim())    error "Samplesheet has a row with an empty 'run' value: ${row}"
            if (!row.fastq_2?.trim())
                error "Sample ${row.sample} run ${row.run} has no fastq_2. This pipeline requires paired-end reads."

            def meta = [
                id       : "${row.sample}_${row.run}".toString(),
                sample   : row.sample.trim(),
                run      : row.run.trim(),
                single_end: false
            ]
            [ meta, [ file(row.fastq_1, checkIfExists: true),
                      file(row.fastq_2, checkIfExists: true) ] ]
        }

    // Optional cohort subset for the scaling benchmarks. Subsetting is applied
    // per *sample*, never per run, so a sample never loses half its data.
    if (subset) {
        ch_keep = ch_reads
            .map { meta, reads -> meta.sample }
            .unique()
            .toSortedList()
            .map { it.take(subset as int) }

        ch_reads = ch_reads
            .combine(ch_keep)
            .filter { meta, reads, keep -> keep.contains(meta.sample) }
            .map { meta, reads, keep -> [meta, reads] }
    }

    // Fail loudly on a duplicated sample+run rather than silently overwriting.
    ch_reads
        .map { meta, reads -> meta.id }
        .toList()
        .map { ids ->
            def dups = ids.countBy { it }.findAll { k, v -> v > 1 }.keySet()
            if (dups) error "Duplicated sample+run combinations in samplesheet: ${dups.join(', ')}"
            ids.size()
        }
        .subscribe { n -> log.info "INPUT_CHECK: ${n} run(s) accepted" }

    emit:
    reads = ch_reads
}
