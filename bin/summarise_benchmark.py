#!/usr/bin/env python3
"""
Summarise a Nextflow execution trace into resource evidence.

Reports, per process: task count, wall time, CPU efficiency, peak memory
against what was requested, and a right-sized recommendation drawn from the
95th percentile of observed use. Right-sizing matters for an allocation case
because a request padded "just in case" is both more expensive and less
credible than one derived from measurement.

    bin/summarise_benchmark.py results/benchmarks/trace-*.txt
    bin/summarise_benchmark.py trace.txt --format markdown -o benchmark.md
    bin/summarise_benchmark.py trace.txt --scale-to 161 --scale-from 8
"""
import argparse
import csv
import glob
import math
import re
import statistics as st
import sys
from collections import defaultdict

# ---------------------------------------------------------------- parsing ---

_DUR = re.compile(r"(\d+(?:\.\d+)?)\s*(ms|s|m|h|d)")
_MEM = re.compile(r"(\d+(?:\.\d+)?)\s*([KMGTP]?B)", re.I)
_MULT_T = {"ms": 1e-3, "s": 1, "m": 60, "h": 3600, "d": 86400}
_MULT_M = {"B": 1, "KB": 2**10, "MB": 2**20, "GB": 2**30, "TB": 2**40, "PB": 2**50}


def dur_s(v):
    """Nextflow durations look like '2h 31m 4s'; sum every component."""
    if not v or v in ("-", "0"):
        return 0.0
    total, found = 0.0, False
    for num, unit in _DUR.findall(str(v)):
        total += float(num) * _MULT_T[unit.lower()]
        found = True
    return total if found else 0.0


def mem_b(v):
    if not v or v in ("-", "0"):
        return 0.0
    m = _MEM.search(str(v))
    if not m:
        return 0.0
    return float(m.group(1)) * _MULT_M.get(m.group(2).upper(), 1)


def pct(v):
    if not v or v in ("-",):
        return 0.0
    try:
        return float(str(v).replace("%", ""))
    except ValueError:
        return 0.0


def p95(xs):
    if not xs:
        return 0.0
    s = sorted(xs)
    return s[min(len(s) - 1, int(math.ceil(0.95 * len(s)) - 1))]


def human_mem(b):
    if b <= 0:
        return "0"
    for unit, size in (("TB", 2**40), ("GB", 2**30), ("MB", 2**20)):
        if b >= size:
            return f"{b / size:.1f} {unit}"
    return f"{b / 1024:.0f} KB"


def human_time(s):
    if s <= 0:
        return "0s"
    if s < 60:
        return f"{s:.0f}s"
    if s < 3600:
        return f"{s / 60:.1f}m"
    if s < 86400:
        return f"{s / 3600:.2f}h"
    return f"{s / 86400:.2f}d"


def load(paths):
    rows = []
    for p in paths:
        with open(p) as fh:
            head = fh.readline()
            delim = "\t" if "\t" in head else ","
            fh.seek(0)
            for r in csv.DictReader(fh, delimiter=delim):
                r["_source"] = p
                rows.append(r)
    return rows


# ------------------------------------------------------------- aggregating --

def summarise(rows, only_ok=True):
    by = defaultdict(list)
    for r in rows:
        if only_ok and r.get("status") != "COMPLETED":
            continue
        by[r.get("process", "?")].append(r)

    out = []
    for proc, ts in sorted(by.items()):
        cpus  = [float(t["cpus"]) for t in ts if t.get("cpus", "").strip().isdigit()]
        rt    = [dur_s(t.get("realtime")) for t in ts]
        pcpu  = [pct(t.get("%cpu")) for t in ts]
        prss  = [mem_b(t.get("peak_rss")) for t in ts]
        reqm  = [mem_b(t.get("memory")) for t in ts]
        rchar = [mem_b(t.get("rchar")) for t in ts]
        wchar = [mem_b(t.get("wchar")) for t in ts]

        ncpu = st.median(cpus) if cpus else 1
        # CPU efficiency: %cpu is reported against a single core, so 100% x
        # ncpu is perfect utilisation of the request.
        eff = (st.mean(pcpu) / (100.0 * ncpu) * 100.0) if (pcpu and ncpu) else 0.0
        req_mem = st.median(reqm) if reqm else 0.0
        peak95  = p95(prss)
        mem_eff = (peak95 / req_mem * 100.0) if req_mem else 0.0

        # Recommendation: p95 peak plus 25% headroom, floored at 1 GB and
        # rounded up to a whole GB.
        rec_mem = max(1.0, math.ceil((peak95 * 1.25) / 2**30)) if peak95 else None
        rec_cpu = ncpu
        if eff and eff < 40 and ncpu > 1:
            rec_cpu = max(1, int(round(ncpu * max(eff, 10) / 100.0 * 1.3)))

        out.append({
            "process": proc,
            "n": len(ts),
            "cpus": ncpu,
            "cpu_eff_pct": eff,
            "realtime_med": st.median(rt) if rt else 0,
            "realtime_max": max(rt) if rt else 0,
            "realtime_sum": sum(rt),
            "cpu_hours": sum(rt) * ncpu / 3600.0,
            "req_mem": req_mem,
            "peak_rss_med": st.median(prss) if prss else 0,
            "peak_rss_p95": peak95,
            "peak_rss_max": max(prss) if prss else 0,
            "mem_eff_pct": mem_eff,
            "rec_cpus": rec_cpu,
            "rec_mem_gb": rec_mem,
            "read_gb": sum(rchar) / 2**30,
            "write_gb": sum(wchar) / 2**30,
        })
    return out


def failures(rows):
    bad = defaultdict(lambda: defaultdict(int))
    for r in rows:
        s = r.get("status")
        if s and s != "COMPLETED":
            bad[r.get("process", "?")][s] += 1
    return bad


# ---------------------------------------------------------------- output ----

def render_text(sums, fails, scale=None):
    L = []
    L.append("=" * 108)
    L.append("RESOURCE SUMMARY BY PROCESS")
    L.append("=" * 108)
    L.append(f"{'process':<32}{'n':>5}{'cpu':>5}{'cpu%':>7}"
             f"{'med time':>10}{'max time':>10}{'req mem':>10}"
             f"{'p95 rss':>10}{'mem%':>7}{'rec':>12}")
    L.append("-" * 108)
    for s in sums:
        rec = f"{int(s[chr(39)+chr(39)] if False else s['rec_cpus']):d}c/{s['rec_mem_gb']}G" if s["rec_mem_gb"] else "-"
        L.append(f"{s['process'][:31]:<32}{s['n']:>5}{s['cpus']:>5.0f}"
                 f"{s['cpu_eff_pct']:>6.0f}%"
                 f"{human_time(s['realtime_med']):>10}"
                 f"{human_time(s['realtime_max']):>10}"
                 f"{human_mem(s['req_mem']):>10}"
                 f"{human_mem(s['peak_rss_p95']):>10}"
                 f"{s['mem_eff_pct']:>6.0f}%{rec:>12}")
    L.append("-" * 108)
    tot_cpuh = sum(s["cpu_hours"] for s in sums)
    tot_task = sum(s["n"] for s in sums)
    L.append(f"{'TOTAL':<32}{tot_task:>5}{'':>5}{'':>7}{'':>10}{'':>10}"
             f"{'':>10}{'':>10}{'':>7}{tot_cpuh:>10.1f}h")
    L.append("")
    L.append(f"Total CPU-hours       : {tot_cpuh:,.1f}")
    L.append(f"Total tasks           : {tot_task:,}")
    L.append(f"Total read / written  : {sum(s['read_gb'] for s in sums):,.0f} GB "
             f"/ {sum(s['write_gb'] for s in sums):,.0f} GB")

    if fails:
        L.append("")
        L.append("NON-COMPLETED TASKS")
        for proc, st_counts in sorted(fails.items()):
            for status, n in sorted(st_counts.items()):
                L.append(f"  {proc:<40} {status:<12} {n}")

    if scale:
        factor, frm, to = scale
        L.append("")
        L.append("=" * 108)
        L.append(f"PROJECTION  {frm} -> {to} samples  (x{factor:.2f})")
        L.append("=" * 108)
        L.append(f"{'process':<32}{'proj CPU-h':>14}{'proj tasks':>14}")
        L.append("-" * 108)
        for s in sums:
            L.append(f"{s['process'][:31]:<32}"
                     f"{s['cpu_hours'] * factor:>14,.1f}"
                     f"{s['n'] * factor:>14,.0f}")
        L.append("-" * 108)
        L.append(f"{'PROJECTED TOTAL':<32}{tot_cpuh * factor:>14,.1f}")
        L.append("")
        L.append("Scaling is linear in sample count, which holds for the per-sample")
        L.append("stages. Cohort-wide steps (GenomicsDBImport, GenotypeGVCFs) scale")
        L.append("super-linearly in memory and are better measured than extrapolated.")
    return "\n".join(L)


def render_markdown(sums, fails, scale=None):
    L = ["# Benchmark summary", "", "## Resource use by process", "",
         "| Process | Tasks | CPUs | CPU eff | Median time | Max time | "
         "Requested mem | p95 peak RSS | Mem eff | Recommended |",
         "|---|--:|--:|--:|--:|--:|--:|--:|--:|---|"]
    for s in sums:
        rec = f"`{int(s['rec_cpus'])} cpus / {int(s['rec_mem_gb'])} GB`" if s["rec_mem_gb"] else "-"
        L.append(f"| `{s['process']}` | {s['n']} | {s['cpus']:.0f} | "
                 f"{s['cpu_eff_pct']:.0f}% | {human_time(s['realtime_med'])} | "
                 f"{human_time(s['realtime_max'])} | {human_mem(s['req_mem'])} | "
                 f"{human_mem(s['peak_rss_p95'])} | {s['mem_eff_pct']:.0f}% | {rec} |")
    tot = sum(s["cpu_hours"] for s in sums)
    L += ["", f"**Total CPU-hours:** {tot:,.1f}  ",
          f"**Total tasks:** {sum(s['n'] for s in sums):,}  ",
          f"**I/O:** {sum(s['read_gb'] for s in sums):,.0f} GB read, "
          f"{sum(s['write_gb'] for s in sums):,.0f} GB written"]
    if fails:
        L += ["", "## Non-completed tasks", "", "| Process | Status | Count |",
              "|---|---|--:|"]
        for proc, sc in sorted(fails.items()):
            for status, n in sorted(sc.items()):
                L.append(f"| `{proc}` | {status} | {n} |")
    if scale:
        factor, frm, to = scale
        L += ["", f"## Projection: {frm} to {to} samples (x{factor:.2f})", "",
              "| Process | Projected CPU-hours | Projected tasks |", "|---|--:|--:|"]
        for s in sums:
            L.append(f"| `{s['process']}` | {s['cpu_hours'] * factor:,.1f} | "
                     f"{s['n'] * factor:,.0f} |")
        L += ["", f"**Projected total CPU-hours:** {tot * factor:,.1f}"]
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("trace", nargs="+", help="trace file(s); globs allowed")
    ap.add_argument("-o", "--output")
    ap.add_argument("--format", choices=["text", "markdown", "csv"], default="text")
    ap.add_argument("--scale-from", type=float, help="sample count this trace represents")
    ap.add_argument("--scale-to", type=float, help="sample count to project to")
    ap.add_argument("--include-failed", action="store_true",
                    help="include non-completed tasks in the statistics")
    args = ap.parse_args()

    paths = []
    for p in args.trace:
        hits = sorted(glob.glob(p))
        paths.extend(hits if hits else [p])
    if not paths:
        sys.exit("ERROR: no trace files matched")

    rows = load(paths)
    if not rows:
        sys.exit("ERROR: trace files contained no task rows")
    print(f"# {len(rows)} task rows from {len(paths)} trace file(s)", file=sys.stderr)

    sums = summarise(rows, only_ok=not args.include_failed)
    fails = failures(rows)

    scale = None
    if args.scale_to and args.scale_from:
        scale = (args.scale_to / args.scale_from, int(args.scale_from), int(args.scale_to))

    if args.format == "csv":
        import io
        buf = io.StringIO()
        w = csv.DictWriter(buf, fieldnames=list(sums[0].keys()))
        w.writeheader()
        w.writerows(sums)
        text = buf.getvalue()
    elif args.format == "markdown":
        text = render_markdown(sums, fails, scale)
    else:
        text = render_text(sums, fails, scale)

    if args.output:
        open(args.output, "w").write(text + "\n")
        print(f"wrote {args.output}", file=sys.stderr)
    else:
        print(text)


if __name__ == "__main__":
    main()
