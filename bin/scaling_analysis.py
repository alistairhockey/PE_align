#!/usr/bin/env python3
"""
Compare resource use across cohort sizes to measure how each stage scales.

The question an allocation reviewer will ask is not "what did it cost" but
"what will it cost at full size". Those differ by stage:

  per-sample stages   FASTP, BWA_MEM, MARKDUPLICATES, HAPLOTYPECALLER.
                      Work is proportional to samples, so cost scales ~linearly
                      (exponent ~1.0) and extrapolation is safe.

  cohort-wide stages  GENOMICSDBIMPORT, GENOTYPEGVCFS. Every sample's data for
                      one interval is held at once, so both time and especially
                      memory can grow faster than linearly. Extrapolating these
                      from a small cohort understates them, which is the single
                      biggest risk in a resource request.

Rather than assume linearity, this fits cost = a * N^b per stage across the
runs you supply and reports the exponent b. Projections then use the measured
exponent.

    bin/scaling_analysis.py \\
        --run 8:results/bench8/benchmarks/trace-bench8-*.txt \\
        --run 24:results/bench24/benchmarks/trace-bench24-*.txt \\
        --run 161:results/bench_238/benchmarks/trace-bench_238-*.txt \\
        --project-to 161

Each --run is <sample_count>:<trace file>. Sample count is the number of
BIOLOGICAL samples, not sequencing runs -- cohort-wide memory scales with
samples, since that is what GenomicsDBImport holds concurrently.
"""
import argparse
import csv
import glob
import math
import re
import statistics as st
import sys
from collections import defaultdict

STAGES = [
    ("Read download and QC (fastp)",
     ["SRA_FETCH", "SRA_RESOLVE", "FASTP", "FASTQC"], "per-sample"),
    ("Alignment (bwa) and duplicate marking",
     ["BWA_INDEX", "BWA_MEM", "SAMTOOLS_MERGE", "MARKDUPLICATES",
      "SAMTOOLS_INDEX", "SAMTOOLS_STATS", "NORMALISE_FASTA", "SAMTOOLS_FAIDX",
      "CREATESEQUENCEDICTIONARY", "BUILD_INTERVALS"], "per-sample"),
    ("Per-sample variant calling (scattered)",
     ["HAPLOTYPECALLER"], "per-sample"),
    ("Joint genotyping and filtering",
     ["GENOMICSDBIMPORT", "GENOTYPEGVCFS", "BCFTOOLS_CONCAT",
      "VARIANTFILTRATION", "MULTIQC"], "cohort-wide"),
    ("Sliding-window and locus PCA",
     ["VCF_PREP_PCA", "VCF_SPLIT_CHROM", "VCF_SUBSET_LOCI"], "cohort-wide"),
]


def dur(v):
    if not v or v == "-":
        return 0.0
    return sum(float(n) * {"ms": 1e-3, "s": 1, "m": 60, "h": 3600, "d": 86400}[u]
               for n, u in re.findall(r"(\d+(?:\.\d+)?)\s*(ms|s|m|h|d)", str(v)))


def mem(v):
    if not v or v == "-":
        return 0.0
    m = re.match(r"^([\d.]+)\s*([KMGTP]?B)", str(v).strip(), re.I)
    if not m:
        return 0.0
    return float(m.group(1)) * {"B": 1, "KB": 2**10, "MB": 2**20,
                                "GB": 2**30, "TB": 2**40}[m.group(2).upper()]


def stage_of(proc):
    for name, keys, kind in STAGES:
        if any(k in proc for k in keys):
            return name, kind
    return None, None


def load(path):
    """Per-stage totals from one trace. CACHED tasks count: their metrics are real."""
    hits = sorted(glob.glob(path))
    if not hits:
        sys.exit(f"no trace matched: {path}")
    if len(hits) > 1:
        print(f"# {path} matched {len(hits)} files, using newest: {hits[-1]}",
              file=sys.stderr)
    f = hits[-1]

    agg = defaultdict(lambda: dict(tasks=0, corehours=0.0, wall=0.0,
                                   peak=0.0, rss=[]))
    with open(f) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            if (r.get("status") or "") not in ("COMPLETED", "CACHED"):
                continue
            proc = r.get("process") or ""
            name, _ = stage_of(proc.split(":")[-1])
            if not name:
                continue
            try:
                cpus = int(r.get("cpus") or 1)
            except ValueError:
                cpus = 1
            rt = dur(r.get("realtime"))
            pk = mem(r.get("peak_rss"))
            a = agg[name]
            a["tasks"] += 1
            a["wall"] += rt
            a["corehours"] += rt * cpus / 3600.0
            a["peak"] = max(a["peak"], pk)
            if pk:
                a["rss"].append(pk)
    return f, agg


def fit_exponent(xs, ys):
    """Least-squares slope of log(y) vs log(x): the exponent b in y = a*N^b."""
    pts = [(x, y) for x, y in zip(xs, ys) if x > 0 and y > 0]
    if len(pts) < 2:
        return None, None
    lx = [math.log(x) for x, _ in pts]
    ly = [math.log(y) for _, y in pts]
    n = len(pts)
    mx, my = sum(lx) / n, sum(ly) / n
    den = sum((x - mx) ** 2 for x in lx)
    if den == 0:
        return None, None
    b = sum((x - mx) * (y - my) for x, y in zip(lx, ly)) / den
    a = math.exp(my - b * mx)
    return b, a


def verdict(b, kind):
    if b is None:
        return "insufficient data"
    if b < 0.8:
        return "sub-linear"
    if b <= 1.2:
        return "linear"
    if b <= 1.6:
        return "super-linear"
    return "STEEP - measure, do not extrapolate"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True,
                    metavar="N:TRACE", help="samples:trace-file (repeatable)")
    ap.add_argument("--project-to", type=float,
                    help="project each stage to this sample count")
    ap.add_argument("-o", "--output")
    args = ap.parse_args()

    runs = []
    for spec in args.run:
        if ":" not in spec:
            sys.exit(f"--run needs <samples>:<trace>, got {spec!r}")
        n, path = spec.split(":", 1)
        f, agg = load(path)
        runs.append((float(n), f, agg))
    runs.sort(key=lambda r: r[0])

    L = []
    L.append("=" * 100)
    L.append("SCALING BY STAGE")
    L.append("=" * 100)
    for n, f, _ in runs:
        L.append(f"  {n:>6.0f} samples : {f}")
    L.append("")

    sizes = [n for n, _, _ in runs]
    for name, _, kind in STAGES:
        ch = [a.get(name, {}).get("corehours", 0.0) for _, _, a in runs]
        tk = [a.get(name, {}).get("tasks", 0) for _, _, a in runs]
        pk = [a.get(name, {}).get("peak", 0.0) for _, _, a in runs]
        if not any(ch):
            L.append(f"{name}\n    no tasks recorded in any run")
            L.append("")
            continue

        b_ch, a_ch = fit_exponent(sizes, ch)
        b_pk, _ = fit_exponent(sizes, pk)
        L.append(f"{name}   [{kind}]")
        L.append(f"    {'samples':>9} {'tasks':>8} {'core-h':>10} "
                 f"{'core-h/sample':>15} {'peak RSS':>11}")
        for i, n in enumerate(sizes):
            per = ch[i] / n if n else 0
            L.append(f"    {n:>9.0f} {tk[i]:>8} {ch[i]:>10.1f} {per:>15.2f} "
                     f"{pk[i]/2**30:>10.1f}G")
        if b_ch is not None:
            L.append(f"    cost exponent  b = {b_ch:5.2f}   ({verdict(b_ch, kind)})")
        if b_pk is not None:
            L.append(f"    memory exponent b = {b_pk:5.2f}   ({verdict(b_pk, kind)})")
        if args.project_to and b_ch is not None:
            proj = a_ch * (args.project_to ** b_ch)
            lin = ch[-1] / sizes[-1] * args.project_to
            L.append(f"    projected to {args.project_to:.0f} samples: "
                     f"{proj:,.0f} core-h  (linear assumption would give {lin:,.0f})")
        L.append("")

    if len(sizes) < 3:
        L.append("NOTE: with fewer than three cohort sizes an exponent is")
        L.append("interpolation, not a fit. Treat it as indicative and add a")
        L.append("third point before quoting it.")
    L.append("=" * 100)

    text = "\n".join(L)
    if args.output:
        open(args.output, "w").write(text + "\n")
        print(f"wrote {args.output}", file=sys.stderr)
    else:
        print(text)


if __name__ == "__main__":
    main()
