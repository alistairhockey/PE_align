/*
 * Normalise a reference FASTA before anything indexes it.
 *
 * Assemblies arrive in whatever state they were written in. Two things quietly
 * break a GATK stack:
 *
 *   CRLF line endings   bwa (kseq) and samtools strip the CR, but not every
 *                       reader does, and a CR treated as sequence shifts every
 *                       downstream coordinate. PBA_HatTrick.fasta ships with
 *                       CRLF throughout.
 *   Header whitespace   trailing spaces become part of the sequence name in
 *                       some readers and not others, so the .dict and the BAM
 *                       header disagree.
 *
 * The cleanup is unconditional. An earlier version tried to detect CRLF first
 * and skip the work when the file looked clean; the detection silently failed
 * inside the container and reported a CRLF file as clean, which is a worse
 * outcome than not checking at all. Stripping unconditionally costs one pass
 * and cannot report a false negative.
 *
 * The result is then verified: any CR surviving in the output is a hard error.
 */
process NORMALISE_FASTA {
    tag        "${fasta.name}"
    label      'process_low'
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode

    conda      "conda-forge::gawk=5.3.0"
    container  "quay.io/biocontainers/gawk:5.3.0"

    input:
    path fasta

    output:
    path "normalised/*"   , emit: fasta
    path "fasta_check.txt", emit: report

    script:
    def out = fasta.name.replaceAll(/\.gz\$/, '')
    """
    set -euo pipefail
    mkdir -p normalised

    if [ "${fasta.extension}" = "gz" ]; then
        gzip -cd ${fasta} > tmp.fa
    else
        cat ${fasta} > tmp.fa
    fi
    bytes_in=\$(wc -c < tmp.fa)

    # Strip CR everywhere, then trailing whitespace on header lines only.
    awk '{ gsub(/\\r/, ""); if (/^>/) sub(/[ \\t]+\$/, ""); print }' \\
        tmp.fa > normalised/${out}
    rm -f tmp.fa

    bytes_out=\$(wc -c < normalised/${out})

    # Verify: nothing may survive that we claimed to remove.
    if awk '/\\r/ { found = 1; exit } END { exit !found }' normalised/${out}; then
        echo "ERROR: carriage returns still present after normalisation" >&2
        exit 1
    fi

    # Base count excludes all whitespace, so it is directly comparable to the
    # sum of sequence lengths in the .fai built from this file. A mismatch
    # between the two is the signal that normalisation went wrong.
    bases=\$(grep -v '^>' normalised/${out} | tr -d ' \\t\\n\\r' | wc -c)
    seqs=\$(grep -c '^>' normalised/${out})

    {
      echo "input:          ${fasta.name}"
      echo "output:         ${out}"
      echo "bytes_in:       \$bytes_in"
      echo "bytes_out:      \$bytes_out"
      echo "bytes_removed:  \$(( bytes_in - bytes_out ))"
      echo "crlf_stripped:  \$( [ \$bytes_in -gt \$bytes_out ] && echo yes || echo "none present" )"
      echo "cr_remaining:   0 (verified)"
      echo "sequences:      \$seqs"
      echo "total_bases:    \$bases"
    } > fasta_check.txt

    cat fasta_check.txt
    """
}
