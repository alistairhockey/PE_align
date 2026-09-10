/*
 * Normalise a reference FASTA before anything indexes it.
 *
 * Assemblies arrive in whatever state they were written in. The two things
 * that quietly break downstream tools:
 *
 *   CRLF line endings  - bwa (kseq) and samtools strip the CR, but not every
 *                        tool in a GATK stack does, and a stray CR inside a
 *                        sequence line shifts coordinates. PBA_HatTrick.fasta
 *                        ships with CRLF.
 *   Trailing whitespace in headers - becomes part of the sequence name in some
 *                        readers and not others, so the .dict and the BAM
 *                        header disagree.
 *
 * If the file is already clean this emits a symlink and costs nothing.
 */
process NORMALISE_FASTA {
    tag        "${fasta.name}"
    label      'process_low'
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode

    conda      "conda-forge::sed"
    container  "quay.io/biocontainers/gawk:5.3.0"

    input:
    path fasta

    output:
    path "normalised/*", emit: fasta
    path "fasta_check.txt", emit: report

    script:
    def out = fasta.name.replaceAll(/\.gz$/, '')
    """
    mkdir -p normalised

    if [ "${fasta.extension}" = "gz" ]; then
        gzip -cd ${fasta} > tmp.fa
    else
        cp -P ${fasta} tmp.fa 2>/dev/null || cat ${fasta} > tmp.fa
    fi

    # Does it actually need fixing? Only inspect the head -- a 700 MB scan to
    # answer a yes/no question is wasted work.
    if head -c 2000000 tmp.fa | grep -qU \$'\\r'; then
        crlf=yes
    else
        crlf=no
    fi

    if [ "\$crlf" = "yes" ]; then
        # Strip CR everywhere, and any trailing whitespace on header lines.
        sed -e 's/\\r\$//' -e '/^>/s/[[:space:]]*\$//' tmp.fa > normalised/${out}
        rm -f tmp.fa
    else
        mv tmp.fa normalised/${out}
    fi

    {
      echo "input:        ${fasta.name}"
      echo "output:       ${out}"
      echo "crlf_found:   \$crlf"
      echo "crlf_fixed:   \$( [ "\$crlf" = "yes" ] && echo yes || echo "not needed" )"
      echo "sequences:    \$(grep -c '^>' normalised/${out})"
      echo "total_bases:  \$(grep -v '^>' normalised/${out} | tr -d '\\n' | wc -c)"
    } > fasta_check.txt

    cat fasta_check.txt
    """
}
