/*
 * Derive scatter intervals from the reference index.
 *
 * Replaces the original pipeline's `Channel.fromPath(intervals) | splitText()`,
 * which passed interval names to GATK with their trailing newline still
 * attached and silently included unplaced scaffolds.
 */
process BUILD_INTERVALS {
    tag        "${fai}"
    label      'process_single'
    publishDir "${params.outdir}/reference/intervals", mode: params.publish_dir_mode

    conda      "conda-forge::gawk=5.3.0"
    container  "quay.io/biocontainers/gawk:5.3.0"

    input:
    path fai
    val  min_length
    val  chr_regex

    output:
    path "*.intervals" , emit: intervals
    path "intervals.tsv"  , emit: table

    script:
    def re = chr_regex ?: '.'
    """
    # Plain-text GATK interval files ('.intervals'), one region per line as
    # <contig>:<start>-<end>.
    #
    # NOT Picard '.interval_list' format. That format is parsed by
    # IntervalListCodec, which requires an @SQ sequence-dictionary header at
    # the top of the file to resolve contig names. A headerless file does not
    # error -- every interval is silently dropped as "unknown reference",
    # GATK reports "Processing 0 bp from intervals", and HaplotypeCaller then
    # writes an empty gVCF and exits 0. The plain-text form carries no
    # dictionary requirement and cannot fail that way.
    awk -F'\t' -v minlen=${min_length} -v re='${re}' '
        \$2 >= minlen && \$1 ~ re {
            name = \$1
            gsub(/[^A-Za-z0-9._-]/, "_", name)     # safe for use as a filename
            print \$1":1-"\$2 > (name ".intervals")
            print \$1"\t"\$2 >> "intervals.tsv"
        }
    ' ${fai}

    if [ ! -s intervals.tsv ]; then
        echo "ERROR: no reference sequences matched regex '${re}' with length >= ${min_length}" >&2
        echo "Check --intervals_min_length and --chr_regex for this genome." >&2
        exit 1
    fi

    # Every emitted file must contain exactly one usable region.
    for f in *.intervals; do
        if ! grep -qE '^[^[:space:]]+:[0-9]+-[0-9]+\$' "\$f"; then
            echo "ERROR: \$f is not a valid interval file:" >&2
            cat "\$f" >&2
            exit 1
        fi
    done
    """
}
