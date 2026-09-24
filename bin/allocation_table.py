#!/usr/bin/env python3
"""
Build the stage-level resource table for an HPC allocation application.

Groups the pipeline's processes into the stages a reviewer thinks in, and keeps
the arithmetic internally consistent: CPU h = Jobs x Cores/job x Wall-time,
which is the first thing anyone checks.

    bin/allocation_table.py --since 2026-09-20
    bin/allocation_table.py --since 2026-09-20 --scale-from 24 --scale-to 161
    bin/allocation_table.py --since 2026-09-20 --format markdown -o table.md

Wall-time is the MEAN per job, so the product above holds. Memory is the
observed peak RSS plus 25% headroom, rounded up -- the figure worth requesting,
not the figure that was requested.
"""
import argparse
import math
import os
import statistics as st
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarise_sacct import collect                      # noqa: E402

# Ordered; each stage lists the process-name substrings that belong to it.
STAGES = [
    ("Read download and QC (fastp)",
     ["SRA_FETCH", "SRA_RESOLVE", "FASTP", "FASTQC"]),
    ("Alignment (bwa) and duplicate marking",
     ["BWA_INDEX", "BWA_MEM", "SAMTOOLS_MERGE", "MARKDUPLICATES",
      "SAMTOOLS_INDEX", "SAMTOOLS_STATS", "NORMALISE_FASTA", "SAMTOOLS_FAIDX",
      "CREATESEQUENCEDICTIONARY", "BUILD_INTERVALS"]),
    ("Per-sample variant calling (scattered)",
     ["HAPLOTYPECALLER"]),
    ("Joint genotyping and filtering",
     ["GENOMICSDBIMPORT", "GENOTYPEGVCFS", "BCFTOOLS_CONCAT",
      "VARIANTFILTRATION", "MULTIQC"]),
    ("Sliding-window and locus PCA",
     ["VCF_PREP_PCA", "VCF_SPLIT_CHROM", "VCF_SUBSET_LOCI"]),
    ("Haplotype block characterisation", []),
    ("Breakpoint mapping against new assemblies", []),
]


def stage_of(proc):
    for name, keys in STAGES:
        if any(k in proc for k in keys):
            return name
    return None


def dedupe(rows):
    """
    Collapse retries: one record per logical task, keyed on process + tag.

    sacct holds a row per SLURM job, so a task that ran in several pipeline
    invocations appears several times and inflates both job counts and
    core-hours.

    The record kept is the LATEST by start time, not the first and not the
    longest. That matters here: before the interval-format fix, HaplotypeCaller
    resolved zero bases and completed in ~6 seconds, so those runs are all
    COMPLETED with near-zero wall time. Keeping the first record per task
    reported a 0.1-minute median against a true 14.6 minutes. The most recent
    execution is the one that reflects the current pipeline.

    MaxRSS is still taken as the maximum across every attempt, since that is
    the memory the task actually needed at some point.
    """
    best, attempts = {}, defaultdict(int)
    for r in rows:
        key = (r["process"], r["tag"])
        attempts[key] += 1
        cur = best.get(key)
        if cur is None:
            best[key] = dict(r)
            continue
        peak = max(cur["rss"], r["rss"])
        # Prefer a COMPLETED record; among those, the most recent start.
        take = ((r["state"] == "COMPLETED" and cur["state"] != "COMPLETED") or
                (r["state"] == "COMPLETED") == (cur["state"] == "COMPLETED")
                and (r.get("start") or "") > (cur.get("start") or ""))
        if take:
            best[key] = dict(r)
        best[key]["rss"] = peak
    for key, rec in best.items():
        rec["attempts"] = attempts[key]
    return list(best.values())


def build(rows, scale=1.0):
    by = defaultdict(list)
    unmapped = defaultdict(int)
    for r in rows:
        if r["state"] != "COMPLETED":
            continue
        s = stage_of(r["process"])
        if s is None:
            unmapped[r["process"]] += 1
            continue
        by[s].append(r)

    out = []
    for name, _ in STAGES:
        ts = by.get(name, [])
        if not ts:
            out.append(dict(stage=name, jobs=None))
            continue
        jobs = len(ts) * scale
        # Weight cores by wall time: a stage mixing 16-core alignment with
        # 1-core indexing should report the cores that actually cost something.
        wall = sum(t["elapsed"] for t in ts) / 3600.0
        cores = (sum(t["cpus"] * t["elapsed"] for t in ts) /
                 sum(t["elapsed"] for t in ts)) if wall else 1
        mean_wall = (wall / len(ts))
        rss = [t["rss"] for t in ts if t["rss"] > 0]
        mem = math.ceil(max(rss) * 1.25 / 2**30) if rss else None
        retried = sum(1 for t in ts if t.get("attempts", 1) > 1)
        out.append(dict(stage=name, jobs=jobs, cores=cores,
                        wall=mean_wall,
                        mem=mem,
                        cpuh=jobs * cores * mean_wall,
                        retried=retried, ntasks=len(ts)))
    return out, unmapped


def fmt(v, spec, dash="—"):
    return dash if v is None else format(v, spec)


def render(rows, markdown, scale, frm, to):
    tot_jobs = sum(r["jobs"] for r in rows if r.get("jobs"))
    tot_cpuh = sum(r["cpuh"] for r in rows if r.get("jobs"))

    if markdown:
        L = ["| Stage | Jobs | Cores/job | Wall-time (h) | Memory (GB) | CPU h |",
             "|---|--:|--:|--:|--:|--:|"]
        for r in rows:
            if not r.get("jobs"):
                L.append(f"| {r['stage']} | — | — | — | — | — |")
            else:
                L.append(f"| {r['stage']} | {r['jobs']:,.0f} | {r['cores']:.0f} | "
                         f"{r['wall']:.2f} | {fmt(r['mem'],'d')} | {r['cpuh']:,.0f} |")
        L.append(f"| **Total** | **{tot_jobs:,.0f}** | | | | **{tot_cpuh:,.0f}** |")
        return "\n".join(L)

    w = 44
    L = []
    hdr = (f"{'Stage':<{w}}{'Jobs':>8}{'Cores/job':>11}"
           f"{'Wall-time (h)':>15}{'Memory (GB)':>13}{'CPU h':>10}")
    L.append("=" * len(hdr)); L.append(hdr); L.append("-" * len(hdr))
    for r in rows:
        if not r.get("jobs"):
            L.append(f"{r['stage']:<{w}}{'—':>8}{'—':>11}{'—':>15}{'—':>13}{'—':>10}")
        else:
            L.append(f"{r['stage']:<{w}}{r['jobs']:>8,.0f}{r['cores']:>11.0f}"
                     f"{r['wall']:>15.2f}{fmt(r['mem'],'d'):>13}{r['cpuh']:>10,.0f}")
    L.append("-" * len(hdr))
    L.append(f"{'Total':<{w}}{tot_jobs:>8,.0f}{'':>11}{'':>15}{'':>13}{tot_cpuh:>10,.0f}")
    L.append("=" * len(hdr))
    if scale != 1.0:
        L.append(f"Scaled x{scale:.2f} from {frm} to {to} samples.")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--since", required=True)
    ap.add_argument("--user", default=None)
    ap.add_argument("--scale-from", type=float)
    ap.add_argument("--scale-to", type=float)
    ap.add_argument("--format", choices=["text", "markdown", "csv"], default="text")
    ap.add_argument("-o", "--output")
    args = ap.parse_args()

    user = args.user or os.environ.get("USER")
    rows = collect(args.since, user, r"^nf-")
    if not rows:
        sys.exit(f"no nf- jobs in accounting since {args.since}")
    raw = len(rows)
    rows = dedupe(rows)
    print(f"# {raw} sacct records -> {len(rows)} distinct tasks after collapsing retries",
          file=sys.stderr)

    scale, frm, to = 1.0, None, None
    if args.scale_from and args.scale_to:
        scale = args.scale_to / args.scale_from
        frm, to = int(args.scale_from), int(args.scale_to)

    table, unmapped = build(rows, scale)
    if args.format == "csv":
        import io, csv as _csv
        buf = io.StringIO(); w = _csv.writer(buf)
        w.writerow(["Stage","Jobs","Cores/job","Wall-time (h)","Memory (GB)","CPU h"])
        for r in table:
            if not r.get("jobs"):
                w.writerow([r["stage"], "", "", "", "", ""])
            else:
                w.writerow([r["stage"], f"{r['jobs']:.0f}", f"{r['cores']:.0f}",
                            f"{r['wall']:.2f}", r["mem"] if r["mem"] else "",
                            f"{r['cpuh']:.0f}"])
        w.writerow(["Total",
                    f"{sum(r['jobs'] for r in table if r.get('jobs')):.0f}", "", "", "",
                    f"{sum(r['cpuh'] for r in table if r.get('jobs')):.0f}"])
        text = buf.getvalue().rstrip("\n")
    else:
        text = render(table, args.format == "markdown", scale, frm, to)
    if args.output:
        open(args.output, "w").write(text + "\n")
        print(f"wrote {args.output}", file=sys.stderr)
    else:
        print(text)
    if unmapped:
        print("\nWARNING: processes not mapped to any stage: "
              + ", ".join(sorted(unmapped)), file=sys.stderr)


if __name__ == "__main__":
    main()
