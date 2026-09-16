/*
 * Resolve an SRA run table into read pairs, downloading only what is missing.
 *
 * Files already present in --reads_dir are adopted as-is, so a partially
 * completed download (or one done outside the pipeline with bin/fetch_reads.sh)
 * is picked up without re-transfer. Naming is fixed by convention:
 *     <RUN>_<SAMPLE>_1.fastq.gz / _2.fastq.gz
 */

include { SRA_RESOLVE } from '../../modules/local/sra_resolve'
include { SRA_FETCH   } from '../../modules/local/sra_fetch'

workflow FETCH_READS {

    take:
    runtable      // path to NCBI SRA run table
    reads_dir     // path or null: existing FASTQs to reuse
    do_download   // bool
    layout        // val
    assay         // val
    organism      // val
    subset        // int or null

    main:
    SRA_RESOLVE(runtable, layout, assay, organism)

    ch_runs = SRA_RESOLVE.out.manifest
        .splitCsv(header: true, sep: '\t', strip: true)
        .map { row ->
            def meta = [
                id        : "${row.sample}_${row.run}".toString(),
                sample    : row.sample,
                run       : row.run,
                population: row.population,
                single_end: false
            ]
            [ meta, row.url_1, row.md5_1, row.url_2, row.md5_2 ]
        }

    // Subset by sample for the scaling benchmarks, never by run: a sample must
    // keep all of its runs or its coverage is wrong.
    if (subset) {
        // Subset by SAMPLE, never by run: a sample must keep every one of its
        // runs or its coverage is silently wrong.
        //
        // Done as one collect/flatMap rather than combine() against a list
        // channel. combine() spreads a List emission across separate tuple
        // elements, so the filter closure was handed one argument per sample
        // name instead of a single list -- "Invalid method invocation `call`
        // with arguments: [...] on _closure6 type".
        ch_runs = ch_runs
            .toList()
            .flatMap { rows ->
                def keep = rows.collect { it[0].sample }
                                .unique()
                                .sort()
                                .take(subset as int)
                rows.findAll { keep.contains(it[0].sample) }
            }
    }

    // Adopt anything already on disk.
    def rd = reads_dir ? file(reads_dir) : null
    ch_split = ch_runs.branch { meta, u1, m1, u2, m2 ->
        def f1 = rd ? rd.resolve("${meta.run}_${meta.sample}_1.fastq.gz") : null
        def f2 = rd ? rd.resolve("${meta.run}_${meta.sample}_2.fastq.gz") : null
        present: rd && f1.exists() && f2.exists()
            return [ meta, [ f1, f2 ] ]
        absent : true
            return [ meta, u1, m1, u2, m2 ]
    }

    if (!do_download) {
        ch_split.absent
            .map { meta, u1, m1, u2, m2 -> meta.run }
            .collect()
            .subscribe { runs ->
                if (runs)
                    error "--download_reads false, but ${runs.size()} run(s) are " +
                          "missing from --reads_dir: ${runs.take(5).join(', ')}" +
                          (runs.size() > 5 ? ' ...' : '')
            }
        ch_reads = ch_split.present
    } else {
        SRA_FETCH(ch_split.absent)
        ch_reads = ch_split.present
            .mix( SRA_FETCH.out.reads.map { meta, r1, r2 -> [ meta, [r1, r2] ] } )
    }

    ch_reads.count().subscribe { n -> log.info "FETCH_READS: ${n} run(s) ready" }

    emit:
    reads    = ch_reads
    manifest = SRA_RESOLVE.out.manifest
}
