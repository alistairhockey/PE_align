#!/usr/bin/env python3
"""
Build a --bam_input samplesheet from a directory of per-run BAMs.

Filenames written by this pipeline's BWA_MEM encode both identifiers:

    <sample>_<run>.bam        e.g. Bari1_002_SRR6242348.bam

so sample and run are recovered by splitting on the final run accession.
Verifies each BAM is non-empty and, when samtools is available, that its @RG SM
tag matches the sample parsed from the filename -- everything downstream reads
sample identity from the BAM header, not from this sheet, so a mismatch there
would silently mis-assign genotypes.

    bin/make_bam_samplesheet.py -d <bam dir> -o bams.csv
    bin/make_bam_samplesheet.py -d <bam dir> -o bams.csv --manifest assets/cret_runs.tsv
"""
import argparse
import csv
import glob
import os
import re
import subprocess
import sys

RUN_RE = re.compile(r"^(?P<sample>.+)_(?P<run>(?:SRR|ERR|DRR)\d+)\.bam$")


def rg_sample(bam):
    """Read the @RG SM tag, if samtools is on PATH."""
    try:
        out = subprocess.run(["samtools", "view", "-H", bam],
                             capture_output=True, text=True, timeout=120)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    if out.returncode != 0:
        return None
    for line in out.stdout.splitlines():
        if line.startswith("@RG"):
            for f in line.split("\t"):
                if f.startswith("SM:"):
                    return f[3:]
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-d", "--bam-dir", required=True)
    ap.add_argument("-o", "--output", default="bams.csv")
    ap.add_argument("--manifest", help="cret_runs.tsv, to check for missing runs")
    ap.add_argument("--no-header-check", action="store_true",
                    help="skip the @RG SM verification (faster)")
    args = ap.parse_args()

    bams = sorted(glob.glob(os.path.join(args.bam_dir, "*.bam")))
    if not bams:
        sys.exit(f"no BAMs found in {args.bam_dir}")

    rows, skipped, mismatched, empty = [], [], [], []
    for b in bams:
        name = os.path.basename(b)
        m = RUN_RE.match(name)
        if not m:
            skipped.append(name)
            continue
        if os.path.getsize(b) == 0:
            empty.append(name)
            continue
        sample, run = m.group("sample"), m.group("run")
        if not args.no_header_check:
            sm = rg_sample(b)
            if sm and sm != sample:
                mismatched.append(f"{name}: filename says {sample}, @RG SM says {sm}")
                continue
        rows.append({"sample": sample, "run": run, "bam": os.path.abspath(b)})

    for label, items in (("unparseable filename", skipped),
                         ("zero bytes", empty),
                         ("@RG SM mismatch", mismatched)):
        if items:
            print(f"WARNING: {len(items)} skipped [{label}]:", file=sys.stderr)
            for i in items[:5]:
                print(f"    {i}", file=sys.stderr)

    if not rows:
        sys.exit("ERROR: no usable BAMs")

    with open(args.output, "w", newline="") as fh:
        # lineterminator="\n": csv.writer defaults to \r\n, which puts a
        # trailing CR on the last field of every row. Nextflow's splitCsv
        # strip:true happens to absorb it, but any shell or awk consumer of
        # this file gets a path ending in CR -- the same CRLF class of bug as
        # the reference FASTA.
        w = csv.DictWriter(fh, fieldnames=["sample", "run", "bam"],
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    samples = sorted({r["sample"] for r in rows})
    print(f"wrote {args.output}: {len(samples)} samples, {len(rows)} runs",
          file=sys.stderr)

    if args.manifest:
        want = {r["run"] for r in csv.DictReader(open(args.manifest), delimiter="\t")}
        have = {r["run"] for r in rows}
        missing = sorted(want - have)
        if missing:
            print(f"WARNING: {len(missing)} run(s) in the manifest have no BAM: "
                  f"{', '.join(missing[:8])}{' ...' if len(missing) > 8 else ''}",
                  file=sys.stderr)
        else:
            print(f"all {len(want)} manifest runs have a BAM", file=sys.stderr)


if __name__ == "__main__":
    main()
