#!/usr/bin/env python3
"""
Per-process CPU and memory statistics from SLURM accounting.

A trace-independent view of the same information `summarise_benchmark.py`
reports. Use it when the Nextflow trace is unavailable or incomplete -- for
example when a trace file already existed and `trace.overwrite` was false, in
which case Nextflow disables the trace writer and records nothing.

sacct is also the more authoritative source for memory: MaxRSS comes from the
cgroup that actually killed the task, not from Nextflow's polling, which can
miss a short spike.

    bin/summarise_sacct.py --since 2026-09-22
    bin/summarise_sacct.py --since 2026-09-22 --format markdown -o stats.md

Task jobs are matched by Nextflow's naming convention, nf-<PROCESS>_(<tag>).
"""
import argparse
import re
import statistics as st
import subprocess
import sys
from collections import defaultdict

FIELDS = ("JobID,JobName%200,State,Elapsed,TotalCPU,AllocCPUS,"
          "ReqMem,MaxRSS,MaxVMSize,MaxDiskRead,MaxDiskWrite,Start")


def to_bytes(v):
    v = (v or "").strip()
    if not v or v in ("-", "0"):
        return 0.0
    m = re.match(r"^([\d.]+)\s*([KMGTP]?)", v, re.I)
    if not m:
        return 0.0
    mult = {"": 1, "K": 2**10, "M": 2**20, "G": 2**30, "T": 2**40, "P": 2**50}
    return float(m.group(1)) * mult[m.group(2).upper()]


def to_secs(v):
    """SLURM durations: [DD-]HH:MM:SS[.mmm] or MM:SS.mmm"""
    v = (v or "").strip()
    if not v or v == "-":
        return 0.0
    days = 0
    if "-" in v:
        d, v = v.split("-", 1)
        days = int(d)
    parts = v.split(":")
    try:
        parts = [float(p) for p in parts]
    except ValueError:
        return 0.0
    while len(parts) < 3:
        parts.insert(0, 0.0)
    h, m, s = parts[-3:]
    return days * 86400 + h * 3600 + m * 60 + s


def human_b(b):
    if b <= 0:
        return "-"
    for u, s in (("TB", 2**40), ("GB", 2**30), ("MB", 2**20)):
        if b >= s:
            return f"{b/s:.1f} {u}"
    return f"{b/1024:.0f} KB"


def human_t(s):
    if s <= 0:
        return "-"
    if s < 60:
        return f"{s:.0f}s"
    if s < 3600:
        return f"{s/60:.1f}m"
    if s < 86400:
        return f"{s/3600:.2f}h"
    return f"{s/86400:.2f}d"


def p95(xs):
    if not xs:
        return 0.0
    xs = sorted(xs)
    import math
    return xs[min(len(xs) - 1, int(math.ceil(0.95 * len(xs)) - 1))]


def collect(since, user, pattern):
    out = subprocess.run(
        ["sacct", "-u", user, "-S", since, "-P", "-n", "--format=" + FIELDS],
        capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"sacct failed: {out.stderr.strip()}")

    main, steps = {}, defaultdict(dict)
    for line in out.stdout.splitlines():
        f = line.split("|")
        if len(f) < 12:
            continue
        jid, name, state, el, tcpu, cpus, req, rss, vmem, dr, dw, start = f[:12]
        base = jid.split(".")[0]
        if "." in jid:
            # Step record: keep the largest RSS seen for the job.
            cur = steps[base]
            for k, v in (("rss", rss), ("vmem", vmem), ("dr", dr), ("dw", dw)):
                if to_bytes(v) > to_bytes(cur.get(k, "0")):
                    cur[k] = v
            if to_secs(tcpu) > to_secs(cur.get("tcpu", "0")):
                cur["tcpu"] = tcpu
        else:
            main[base] = dict(name=name, state=state, elapsed=el, tcpu=tcpu,
                              cpus=cpus, req=req, start=start)

    rows = []
    for jid, m in main.items():
        if not re.match(pattern, m["name"]):
            continue
        proc = re.sub(r"^nf-", "", m["name"])
        tagm = re.search(r"_\((.*)\)\s*$", m["name"])
        tag  = tagm.group(1) if tagm else ""
        proc = re.sub(r"_\(.*$", "", proc)
        proc = re.sub(r"_$", "", proc)
        s = steps.get(jid, {})
        rows.append(dict(
            process=proc, tag=tag, state=m["state"],
            elapsed=to_secs(m["elapsed"]),
            tcpu=max(to_secs(m["tcpu"]), to_secs(s.get("tcpu", "0"))),
            cpus=int(m["cpus"] or 1),
            req=to_bytes(m["req"]),
            rss=to_bytes(s.get("rss", "0")),
            vmem=to_bytes(s.get("vmem", "0")),
            read=to_bytes(s.get("dr", "0")),
            write=to_bytes(s.get("dw", "0")),
            start=m.get("start", "")))
    return rows


def summarise(rows, only_ok=True):
    by = defaultdict(list)
    for r in rows:
        if only_ok and r["state"] != "COMPLETED":
            continue
        by[r["process"]].append(r)

    out = []
    for proc, ts in sorted(by.items()):
        el = [t["elapsed"] for t in ts]
        cpus = st.median([t["cpus"] for t in ts]) or 1
        # CPU efficiency: TotalCPU is aggregate core-seconds actually consumed,
        # against (wall time x allocated cores) that were reserved.
        reserved = sum(t["elapsed"] * t["cpus"] for t in ts)
        used = sum(t["tcpu"] for t in ts)
        eff = (used / reserved * 100) if reserved else 0
        rss = [t["rss"] for t in ts if t["rss"] > 0]
        req = st.median([t["req"] for t in ts if t["req"] > 0] or [0])
        pk = p95(rss)
        import math
        rec = max(1, math.ceil(pk * 1.25 / 2**30)) if pk else None
        out.append(dict(process=proc, n=len(ts), cpus=cpus, eff=eff,
                        med=st.median(el), mx=max(el), total=sum(el),
                        corehours=used / 3600,
                        req=req, rss_med=st.median(rss) if rss else 0,
                        rss_p95=pk, rss_max=max(rss) if rss else 0,
                        memeff=(pk / req * 100) if req else 0,
                        rec_gb=rec,
                        read=sum(t["read"] for t in ts),
                        write=sum(t["write"] for t in ts)))
    return out


def render(sums, rows, markdown=False):
    bad = defaultdict(lambda: defaultdict(int))
    for r in rows:
        if r["state"] != "COMPLETED":
            bad[r["process"]][r["state"]] += 1

    if markdown:
        L = ["# CPU and memory usage by process", "",
             "| Process | Tasks | CPUs | CPU eff | Median | Max | Req mem | "
             "p95 RSS | Max RSS | Mem eff | Recommend |",
             "|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|---|"]
        for s in sums:
            L.append(f"| `{s['process']}` | {s['n']} | {s['cpus']:.0f} | "
                     f"{s['eff']:.0f}% | {human_t(s['med'])} | {human_t(s['mx'])} | "
                     f"{human_b(s['req'])} | {human_b(s['rss_p95'])} | "
                     f"{human_b(s['rss_max'])} | {s['memeff']:.0f}% | "
                     f"{s['rec_gb']} GB |" if s['rec_gb'] else "| - |")
        L += ["", f"**Total core-hours:** {sum(s['corehours'] for s in sums):,.1f}"]
        return "\n".join(L)

    L = ["=" * 118,
         "CPU AND MEMORY USAGE BY PROCESS   (source: SLURM accounting)",
         "=" * 118,
         f"{'process':<34}{'n':>4}{'cpu':>5}{'cpu%':>6}{'median':>9}{'max':>9}"
         f"{'req mem':>10}{'p95 RSS':>10}{'max RSS':>10}{'mem%':>6}{'suggest':>10}",
         "-" * 118]
    for s in sums:
        L.append(f"{s['process'][:33]:<34}{s['n']:>4}{s['cpus']:>5.0f}"
                 f"{s['eff']:>5.0f}%{human_t(s['med']):>9}{human_t(s['mx']):>9}"
                 f"{human_b(s['req']):>10}{human_b(s['rss_p95']):>10}"
                 f"{human_b(s['rss_max']):>10}{s['memeff']:>5.0f}%"
                 f"{(str(s['rec_gb'])+' GB') if s['rec_gb'] else '-':>10}")
    L += ["-" * 118,
          f"{'TOTAL':<34}{sum(s['n'] for s in sums):>4}"
          f"{'':>5}{'':>6}{'':>9}{'':>9}{'':>10}{'':>10}{'':>10}{'':>6}"
          f"{sum(s['corehours'] for s in sums):>9,.0f}h",
          "",
          f"Total core-hours consumed : {sum(s['corehours'] for s in sums):,.1f}",
          f"Total wall time in tasks  : {human_t(sum(s['total'] for s in sums))}",
          f"I/O                       : {human_b(sum(s['read'] for s in sums))} read, "
          f"{human_b(sum(s['write'] for s in sums))} written"]
    if bad:
        L += ["", "NON-COMPLETED TASKS"]
        for p, sc in sorted(bad.items()):
            for stt, n in sorted(sc.items()):
                L.append(f"  {p:<44} {stt:<14} {n}")
    L += ["",
          "cpu%    = TotalCPU / (wall x allocated cores). Low means cores idle.",
          "mem%    = p95 MaxRSS / requested. Low means the request is oversized.",
          "suggest = p95 MaxRSS + 25% headroom, rounded up to whole GB.",
          "=" * 118]
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--since", required=True, help="start date, e.g. 2026-09-22")
    ap.add_argument("--user", default=None)
    ap.add_argument("--pattern", default=r"^nf-", help="job-name regex")
    ap.add_argument("--format", choices=["text", "markdown"], default="text")
    ap.add_argument("-o", "--output")
    ap.add_argument("--include-failed", action="store_true")
    args = ap.parse_args()

    import os
    user = args.user or os.environ.get("USER")
    rows = collect(args.since, user, args.pattern)
    if not rows:
        sys.exit(f"no jobs matching {args.pattern!r} since {args.since}")
    sums = summarise(rows, only_ok=not args.include_failed)
    text = render(sums, rows, markdown=(args.format == "markdown"))
    if args.output:
        open(args.output, "w").write(text + "\n")
        print(f"wrote {args.output}", file=sys.stderr)
    else:
        print(text)


if __name__ == "__main__":
    main()
