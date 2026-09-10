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
    container  "biocontainers/gawk:5.3.0"

    input:
    path fai
    val  min_length
    val  chr_regex

    output:
    path "*.interval_list", emit: intervals
    path "intervals.tsv"  , emit: table

    script:
    def re = chr_regex ?: '.'
    """
    # One interval_list per reference sequence: the scatter unit for
    # HaplotypeCaller and GenomicsDBImport.
    awk -F'\\t' -v minlen=${min_length} -v re='${re}' '
        \$2 >= minlen && \$1 ~ re {
            name = \$1
            gsub(/[^A-Za-z0-9._-]/, "_", name)     # safe for use as a filename
            print \$1"\\t1\\t"\$2"\\t+\\t"\$1 > (name ".interval_list")
            print \$1"\\t"\$2 >> "intervals.tsv"
        }
    ' ${fai}

    if [ ! -s intervals.tsv ]; then
        echo "ERROR: no reference sequences matched regex '${re}' with length >= ${min_length}" >&2
        echo "Check --intervals_min_length and the chr_regex for this genome." >&2
        exit 1
    fi
    """
}
