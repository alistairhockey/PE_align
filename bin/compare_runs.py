#!/usr/bin/env python3
"""
Compare CPU and memory use per process between two runs.

Answers the direct question: at a larger cohort, which processes needed more
cores or more memory, and which did not?

    bin/compare_runs.py \\
        --a 24:results/bench24/benchmarks/trace-bench24-*.txt \\
        --b 161:results/bench_238/benchmarks/trace-bench_238-*.txt

Each argument is <label>:<trace file>. CACHED tasks are included -- their
recorded metrics are from a real execution.

Memory is peak RSS (the maximum observed), because that is what determines
whether a task survives. CPU efficiency is mean %cpu against allocated cores;
below ~50% means cores sat idle and the request can come down.
"""
import argparse
import csv
import glob
import re
import statistics as st
import sys
from collections import defaultdict


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


def pct(v):
    try:
        return float(str(v).replace("%", ""))
    except (TypeError, ValueError):
        return 0.0


def load(spec):
    if ":" not in spec:
        sys.exit(f"expected <label>:<trace>, got {spec!r}")
    label, pat = spec.split(":", 1)
    hits = sorted(glob.glob(pat))
    if not hits:
        sys.exit(f"no trace matched: {pat}")
    path = hits[-1]

    agg = defaultdict(lambda: dict(n=0, cpus=[], eff=[], rss=[], wall=[]))
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("status") or "") not in ("COMPLETED", "CACHED"):
                continue
            proc = (r.get("process") or "").split(":")[-1]
            if not proc:
                continue
            try:
                c = int(r.get("cpus") or 1)
            except ValueError:
                c = 1
            a = agg[proc]
            a["n"] += 1
            a["cpus"].append(c)
            p = pct(r.get("%cpu"))
            if p:
                a["eff"].append(p / (100.0 * c) * 100.0)
            m = mem(r.get("peak_rss"))
            if m:
                a["rss"].append(m)
            a["wall"].append(dur(r.get("realtime")))
    return label, path, agg


def g(b):
    return b / 2**30


def arrow(a, b, lo=0.15):
    if a == 0 or b == 0:
        return " "
    ch = (b - a) / a
    if ch > lo:
        return "^"
    if ch < -lo:
        return "v"
    return "="


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", required=True, metavar="LABEL:TRACE")
    ap.add_argument("--b", required=True, metavar="LABEL:TRACE")
    ap.add_argument("--format", choices=["text", "markdown"], default="text")
    ap.add_argument("-o", "--output")
    args = ap.parse_args()

    la, pa, A = load(args.a)
    lb, pb, B = load(args.b)
    procs = sorted(set(A) | set(B))

    L = []
    L.append(f"A = {la:>6} samples   {pa}")
    L.append(f"B = {lb:>6} samples   {pb}")
    L.append("")
    if args.format == "markdown":
        L = [f"**A** = {la} samples (`{pa}`)  ", f"**B** = {lb} samples (`{pb}`)", "",
             f"| Process | Tasks A | Tasks B | Cores | CPU eff A | CPU eff B | "
             f"Peak RSS A | Peak RSS B | RSS change |",
             "|---|--:|--:|--:|--:|--:|--:|--:|--:|"]
        for p in procs:
            a, b = A.get(p), B.get(p)
            ra = max(a["rss"]) if a and a["rss"] else 0
            rb = max(b["rss"]) if b and b["rss"] else 0
            ch = f"{(rb-ra)/ra*100:+.0f}%" if ra and rb else "—"
            L.append(f"| `{p}` | {a['n'] if a else 0} | {b['n'] if b else 0} | "
                     f"{st.median(a['cpus']) if a else (st.median(b['cpus']) if b else 0):.0f} | "
                     f"{st.mean(a['eff']):.0f}%" if a and a['eff'] else "| — ")
        return_text = "\n".join(L)
        if args.output:
            open(args.output, "w").write(return_text + "\n")
        else:
            print(return_text)
        return

    hdr = (f"{'process':<30}{'n A':>5}{'n B':>6}{'cores':>7}"
           f"{'cpu% A':>8}{'cpu% B':>8}{'RSS A':>9}{'RSS B':>9}{'ΔRSS':>9}"
           f"{'wall A':>9}{'wall B':>9}")
    L += ["=" * len(hdr), "CPU AND MEMORY: A vs B", "=" * len(hdr), hdr,
          "-" * len(hdr)]
    for p in procs:
        a, b = A.get(p), B.get(p)
        na, nb = (a["n"] if a else 0), (b["n"] if b else 0)
        cores = st.median(a["cpus"]) if a else (st.median(b["cpus"]) if b else 0)
        ea = st.mean(a["eff"]) if a and a["eff"] else 0
        eb = st.mean(b["eff"]) if b and b["eff"] else 0
        ra = max(a["rss"]) if a and a["rss"] else 0
        rb = max(b["rss"]) if b and b["rss"] else 0
        wa = st.median(a["wall"]) if a and a["wall"] else 0
        wb = st.median(b["wall"]) if b and b["wall"] else 0
        d = f"{(rb-ra)/ra*100:+.0f}%" if ra and rb else "-"
        L.append(f"{p[:29]:<30}{na:>5}{nb:>6}{cores:>7.0f}"
                 f"{ea:>7.0f}%{eb:>7.0f}%{g(ra):>8.1f}G{g(rb):>8.1f}G"
                 f"{d:>9}{wa/60:>8.0f}m{wb/60:>8.0f}m {arrow(ra, rb)}")
    L += ["-" * len(hdr),
          "",
          "ΔRSS is the change in PEAK memory from A to B. A process whose peak",
          "grows with cohort size is cohort-wide (it holds every sample at once);",
          "one that does not is per-sample and its request can stay fixed.",
          "",
          "cpu% is mean utilisation against allocated cores. Below ~50% means",
          "cores sat idle and the core request can come down.",
          "=" * len(hdr)]
    text = "\n".join(L)
    if args.output:
        open(args.output, "w").write(text + "\n")
        print(f"wrote {args.output}", file=sys.stderr)
    else:
        print(text)


if __name__ == "__main__":
    main()
