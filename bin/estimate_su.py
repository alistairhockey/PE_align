#!/usr/bin/env python3
"""
Project measured resource use into Pawsey Setonix service units.

Setonix CPU nodes are dual AMD EPYC 7763: 128 cores and ~230 GB usable RAM,
so ~1.79 GB per core. Two accounting models are reported because they give
materially different numbers and which one applies depends on how the job is
submitted:

  per-core     charged for the cores requested
                   SU = cores x hours
  whole-node   charged for whole nodes; many tasks pack onto one node, so
               this is core-hours divided by a packing efficiency, never a
               whole node per task
                   SU = core-hours / packing_efficiency

A memory-heavy task is charged as if it used the cores that memory implies:
requesting 64 GB on a node with 1.79 GB/core occupies ~36 cores' worth of the
node whether or not they are used. That "memory-implied core count" is usually
what makes a GATK joint-genotyping request more expensive than its core count
suggests, and it is the number worth showing in an allocation case.

    bin/estimate_su.py results/benchmarks/trace-*.txt
    bin/estimate_su.py trace.txt --scale-from 8 --scale-to 161 --format markdown

RATES ARE NOT HARDCODED POLICY. Confirm the current charge rate and partition
accounting for your project with Pawsey before submitting an application;
override with --su-per-core-hour and --node-cores/--node-mem-gb.
"""
import argparse
import math
import sys
from collections import defaultdict

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from summarise_benchmark import (load, summarise, human_time)   # noqa: E402


def project(sums, node_cores, node_mem_gb, su_rate, factor, packing_eff=0.80):
    mem_per_core = node_mem_gb / node_cores
    rows = []
    for s in sums:
        cores = s["cpus"] or 1
        req_gb = s["req_mem"] / 2**30
        # Cores the memory request alone reserves on a node.
        mem_cores = req_gb / mem_per_core if mem_per_core else 0
        eff_cores = max(cores, mem_cores)

        hours = s["realtime_sum"] / 3600.0 * factor
        core_hours = eff_cores * hours
        su_core = core_hours * su_rate
        # Whole-node accounting does NOT mean a node per task: the scheduler
        # packs many tasks onto each node. The realistic penalty is imperfect
        # packing -- stragglers, memory fragmentation, the tail of a scatter.
        su_node = core_hours * su_rate / packing_eff

        rows.append({
            "process": s["process"],
            "tasks": s["n"] * factor,
            "cores": cores,
            "req_gb": req_gb,
            "mem_cores": mem_cores,
            "eff_cores": eff_cores,
            "node_hours_sum": hours,
            "core_hours": core_hours,
            "node_hours": core_hours / node_cores,
            "su_per_core": su_core,
            "su_whole_node": su_node,
            "binding": "memory" if mem_cores > cores else "cores",
        })
    return rows


def render(rows, factor, frm, to, node_cores, node_mem_gb, su_rate, packing_eff=0.80, markdown=False):
    tot_c = sum(r["su_per_core"] for r in rows)
    tot_n = sum(r["su_whole_node"] for r in rows)
    scale_note = (f"scaled x{factor:.2f} from {frm} to {to} samples"
                  if factor != 1 else "measured run, no scaling")

    if markdown:
        L = ["# Setonix service-unit projection", "",
             f"Node model: **{node_cores} cores / {node_mem_gb:.0f} GB** "
             f"({node_mem_gb / node_cores:.2f} GB per core) at "
             f"**{su_rate} SU per core-hour**, packing efficiency "
             f"**{packing_eff:.0%}**  ",
             f"Basis: {scale_note}", "",
             "| Process | Tasks | Cores | Req mem | Mem-implied cores | Binding | "
             "Core-hours | SU (per-core) | SU (whole-node) |",
             "|---|--:|--:|--:|--:|:-:|--:|--:|--:|"]
        for r in sorted(rows, key=lambda x: -x["su_per_core"]):
            L.append(f"| `{r['process']}` | {r['tasks']:,.0f} | {r['cores']:.0f} | "
                     f"{r['req_gb']:.0f} GB | {r['mem_cores']:.1f} | {r['binding']} | "
                     f"{r['core_hours']:,.0f} | "
                     f"{r['su_per_core']:,.0f} | {r['su_whole_node']:,.0f} |")
        L += ["", f"**Total SU (per-core accounting):** {tot_c:,.0f}  ",
              f"**Total SU (whole-node accounting):** {tot_n:,.0f}  ", "",
              "> Confirm the charge rate and partition accounting with Pawsey "
              "before quoting these figures in an application."]
        return "\n".join(L)

    L = ["=" * 104, "SETONIX SERVICE-UNIT PROJECTION", "=" * 104,
         f"node model : {node_cores} cores / {node_mem_gb:.0f} GB "
         f"({node_mem_gb / node_cores:.2f} GB per core)",
         f"charge rate: {su_rate} SU per core-hour",
         f"packing    : {packing_eff:.0%} (whole-node accounting only)",
         f"basis      : {scale_note}", "-" * 104,
         f"{'process':<30}{'tasks':>9}{'cores':>7}{'req mem':>10}"
         f"{'mem-cores':>11}{'bind':>8}{'core-h':>12}{'SU core':>13}{'SU node':>13}"]
    L.append("-" * 104)
    for r in sorted(rows, key=lambda x: -x["su_per_core"]):
        L.append(f"{r['process'][:29]:<30}{r['tasks']:>9,.0f}{r['cores']:>7.0f}"
                 f"{r['req_gb']:>9.0f}G{r['mem_cores']:>11.1f}{r['binding']:>8}"
                 f"{r['core_hours']:>12,.0f}"
                 f"{r['su_per_core']:>13,.0f}{r['su_whole_node']:>13,.0f}")
    tot_ch = sum(r["core_hours"] for r in rows)
    L += ["-" * 104,
          f"{'TOTAL':<30}{'':>9}{'':>7}{'':>10}{'':>11}{'':>8}{tot_ch:>12,.0f}"
          f"{tot_c:>13,.0f}{tot_n:>13,.0f}", "",
          f"Total SU, per-core accounting   : {tot_c:>12,.0f}",
          f"Total SU, whole-node accounting : {tot_n:>12,.0f}", "",
          "Processes marked 'memory' binding reserve more of a node through their",
          "memory request than through their core count. Right-sizing memory on",
          "those is the cheapest way to reduce the allocation.", "",
          "Confirm the current charge rate and partition accounting with Pawsey",
          "before quoting these figures.", "=" * 104]
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("trace", nargs="+")
    ap.add_argument("-o", "--output")
    ap.add_argument("--format", choices=["text", "markdown"], default="text")
    ap.add_argument("--node-cores", type=int, default=128)
    ap.add_argument("--node-mem-gb", type=float, default=230.0)
    ap.add_argument("--su-per-core-hour", type=float, default=1.0)
    ap.add_argument("--scale-from", type=float)
    ap.add_argument("--scale-to", type=float)
    ap.add_argument("--packing-efficiency", type=float, default=0.80,
                    help="fraction of an allocated node actually kept busy (default 0.80)")
    ap.add_argument("--overhead", type=float, default=1.15,
                    help="multiplier for retries, queueing and reruns (default 1.15)")
    args = ap.parse_args()

    import glob
    paths = []
    for p in args.trace:
        hits = sorted(glob.glob(p))
        paths.extend(hits if hits else [p])

    rows = load(paths)
    if not rows:
        sys.exit("ERROR: no task rows found in trace")
    sums = summarise(rows)

    factor = 1.0
    frm = to = None
    if args.scale_from and args.scale_to:
        factor = args.scale_to / args.scale_from
        frm, to = int(args.scale_from), int(args.scale_to)
    factor *= args.overhead

    proj = project(sums, args.node_cores, args.node_mem_gb,
                   args.su_per_core_hour, factor, args.packing_efficiency)
    text = render(proj, factor, frm, to, args.node_cores, args.node_mem_gb,
                  args.su_per_core_hour, args.packing_efficiency,
                  markdown=(args.format == "markdown"))

    if args.output:
        open(args.output, "w").write(text + "\n")
        print(f"wrote {args.output}", file=sys.stderr)
    else:
        print(text)


if __name__ == "__main__":
    main()
