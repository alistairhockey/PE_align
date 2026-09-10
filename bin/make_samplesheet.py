#!/usr/bin/env python3
"""
Build a samplesheet from reads already on disk.

Scans a reads directory against a resolved run manifest and emits a samplesheet
containing only samples whose **every** run is present and complete. That last
part matters: a sample with two runs where only one has downloaded would be
silently half-covered, which corrupts both the benchmark and any PCA.

    bin/make_samplesheet.py -m assets/cret_runs.tsv \\
        -r /group/peg/cicer/cret/reads -o samplesheet.csv

    --check-md5   verify checksums rather than just file sizes (slow but exact)
    --min-runs    require at least this many complete samples, else exit 1
"""
import argparse
import csv
import hashlib
import os
import sys
from collections import defaultdict


def md5_of(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-m", "--manifest", default="assets/cret_runs.tsv")
    ap.add_argument("-r", "--reads-dir", required=True)
    ap.add_argument("-o", "--output", default="samplesheet.csv")
    ap.add_argument("--check-md5", action="store_true")
    ap.add_argument("--min-runs", type=int, default=1)
    ap.add_argument("--max-samples", type=int, default=None,
                    help="keep only the first N complete samples")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.manifest), delimiter="\t"))
    by_sample = defaultdict(list)
    for r in rows:
        by_sample[r["sample"]].append(r)

    complete, partial, missing = {}, [], []
    for sample, runs in by_sample.items():
        ok = []
        for r in runs:
            f1 = os.path.join(args.reads_dir, f"{r['run']}_{r['sample']}_1.fastq.gz")
            f2 = os.path.join(args.reads_dir, f"{r['run']}_{r['sample']}_2.fastq.gz")
            if not (os.path.exists(f1) and os.path.exists(f2)):
                continue
            if os.path.getsize(f1) != int(r["bytes_1"]) or \
               os.path.getsize(f2) != int(r["bytes_2"]):
                continue
            if args.check_md5:
                if md5_of(f1) != r["md5_1"] or md5_of(f2) != r["md5_2"]:
                    print(f"  MD5 mismatch: {r['run']}", file=sys.stderr)
                    continue
            ok.append((r, f1, f2))

        if len(ok) == len(runs):
            complete[sample] = ok
        elif ok:
            partial.append((sample, len(ok), len(runs)))
        else:
            missing.append(sample)

    print(f"samples complete : {len(complete)}", file=sys.stderr)
    print(f"samples partial  : {len(partial)}   (excluded -- would be "
          f"under-covered)", file=sys.stderr)
    print(f"samples missing  : {len(missing)}", file=sys.stderr)
    for s, got, want in partial[:8]:
        print(f"    {s}: {got}/{want} runs", file=sys.stderr)

    keep = sorted(complete)
    if args.max_samples:
        keep = keep[:args.max_samples]

    with open(args.output, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "run", "fastq_1", "fastq_2"])
        n = 0
        for sample in keep:
            for r, f1, f2 in complete[sample]:
                w.writerow([sample, r["run"], os.path.abspath(f1), os.path.abspath(f2)])
                n += 1

    print(f"wrote {args.output}: {len(keep)} sample(s), {n} run(s)", file=sys.stderr)
    if len(keep) < args.min_runs:
        sys.exit(f"ERROR: only {len(keep)} complete sample(s), need {args.min_runs}")


if __name__ == "__main__":
    main()
