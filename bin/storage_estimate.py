#!/usr/bin/env python3
"""
Project the storage an Acacia (object storage) allocation needs.

Measures the artefacts a completed benchmark run actually produced, then scales
them to the full cohort. Splits what belongs in object storage from what is
scratch, because they are requested separately and sized very differently.

    bin/storage_estimate.py --measured-samples 24 --target-samples 161
    bin/storage_estimate.py --measured-samples 24 --target-samples 161 \\
        --format csv -o acacia_storage.csv

Scaling rules, which differ by artefact and are the point of doing this
properly:

  per-sample, linear   BAMs, gVCFs. One set per sample, so size scales with
                       sample count.
  fixed                reference and its indices. Independent of cohort size.
  cohort-wide          the joint-called VCF. Grows with samples, but far more
                       slowly than linearly: variant sites saturate as more
                       samples are added, while each added sample contributes
                       only one more genotype column per site.
"""
import argparse
import csv
import io
import os
import subprocess
import sys

GB = 1024 ** 3


def du_bytes(path):
    if not os.path.exists(path):
        return 0
    out = subprocess.run(["du", "-sb", path], capture_output=True, text=True)
    try:
        return int(out.stdout.split()[0])
    except (IndexError, ValueError):
        return 0


def find_sum(root, pattern):
    out = subprocess.run(
        ["bash", "-c",
         f"find {root} -name '{pattern}' -printf '%s\\n' 2>/dev/null"],
        capture_output=True, text=True)
    vals = [int(x) for x in out.stdout.split() if x.isdigit()]
    return sum(vals), len(vals)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--measured-samples", type=float, required=True)
    ap.add_argument("--target-samples", type=float, required=True)
    ap.add_argument("--reads-dir", default="/group/peg/cicer/cret/reads")
    ap.add_argument("--results-dir", default="results/bench24")
    ap.add_argument("--work-dir", default="/group/sae001/ahockey/PE_align_work")
    ap.add_argument("--save-interval-gvcfs", action="store_true",
                    help="include the per-interval gVCF shards in the estimate")
    ap.add_argument("--format", choices=["text", "csv", "markdown"], default="text")
    ap.add_argument("-o", "--output")
    args = ap.parse_args()

    n0, n1 = args.measured_samples, args.target_samples
    f = n1 / n0

    reads = du_bytes(args.reads_dir)
    bam_m, bam_n = find_sum(args.results_dir, "*.bam")
    if bam_m == 0:
        bam_m, bam_n = find_sum(args.work_dir, "*.bam")
    gv_m, gv_n = find_sum(args.work_dir, "*.g.vcf.gz")
    vcf = du_bytes(os.path.join(args.results_dir, "variants"))
    ref = du_bytes(os.path.join(args.results_dir, "reference"))
    qc = du_bytes(os.path.join(args.results_dir, "qc"))

    per_gvcf = (gv_m / gv_n) if gv_n else 0
    # Merged per-sample gVCF ~= its interval shards summed.
    intervals_per_sample = 8
    gvcf_target = per_gvcf * intervals_per_sample * n1
    shards_target = gvcf_target if args.save_interval_gvcfs else 0

    # BAMs measured per sample, scaled linearly.
    bam_per_sample = bam_m / n0
    bam_target = bam_per_sample * n1

    # Cohort VCF: sites saturate, genotype columns grow linearly. A square-root
    # scaling of the per-sample contribution is a deliberately conservative
    # middle ground -- flag it rather than pretend it is measured.
    vcf_target = vcf * (f ** 0.5)

    rows = [
        # (item, archive?, measured bytes, target bytes, scaling note)
        ("Raw reads (FASTQ, 238 runs)", True, reads, reads,
         "complete, no scaling"),
        ("Per-sample BAM + index", True, bam_m, bam_target,
         f"linear x{f:.2f}"),
        ("Per-sample gVCF", True, per_gvcf * intervals_per_sample * n0,
         gvcf_target, f"linear x{f:.2f}"),
        ("Per-interval gVCF shards", True,
         gv_m if args.save_interval_gvcfs else 0, shards_target,
         "optional" if not args.save_interval_gvcfs else f"linear x{f:.2f}"),
        ("Cohort VCF (raw + filtered)", True, vcf, vcf_target,
         "sub-linear (sites saturate)"),
        ("Reference + indices", True, ref, ref, "fixed"),
        ("QC reports and benchmarks", True, qc, qc * f, f"linear x{f:.2f}"),
    ]
    scratch = [
        ("Nextflow work directory (transient)", False,
         du_bytes(args.work_dir), du_bytes(args.work_dir) * f,
         "delete after publish"),
    ]

    arch_t = sum(r[3] for r in rows)
    scr_t = sum(r[3] for r in scratch)

    if args.format == "csv":
        buf = io.StringIO(); w = csv.writer(buf)
        w.writerow(["Item", "Destination", f"Measured at {n0:.0f} samples (GB)",
                    f"Projected at {n1:.0f} samples (GB)", "Scaling"])
        for it, arch, m, t, note in rows + scratch:
            w.writerow([it, "Acacia" if arch else "Scratch",
                        f"{m/GB:.1f}", f"{t/GB:.1f}", note])
        w.writerow(["TOTAL Acacia (object storage)", "Acacia", "",
                    f"{arch_t/GB:.1f}", ""])
        w.writerow(["TOTAL scratch (peak, transient)", "Scratch", "",
                    f"{scr_t/GB:.1f}", ""])
        text = buf.getvalue().rstrip("\n")
    elif args.format == "markdown":
        L = [f"| Item | Destination | Measured @ {n0:.0f} (GB) | "
             f"Projected @ {n1:.0f} (GB) | Scaling |",
             "|---|---|--:|--:|---|"]
        for it, arch, m, t, note in rows + scratch:
            L.append(f"| {it} | {'Acacia' if arch else 'Scratch'} | "
                     f"{m/GB:,.1f} | {t/GB:,.1f} | {note} |")
        L.append(f"| **Total Acacia** | | | **{arch_t/GB:,.0f}** | |")
        L.append(f"| **Total scratch (peak)** | | | **{scr_t/GB:,.0f}** | |")
        text = "\n".join(L)
    else:
        w1 = 38
        L = ["=" * 96,
             f"STORAGE ESTIMATE   measured at {n0:.0f} samples, projected to {n1:.0f}",
             "=" * 96,
             f"{'Item':<{w1}}{'Dest':>9}{'Measured':>12}{'Projected':>12}  Scaling"]
        L.append("-" * 96)
        for it, arch, m, t, note in rows:
            L.append(f"{it[:w1-1]:<{w1}}{'Acacia':>9}{m/GB:>11.1f}G"
                     f"{t/GB:>11.1f}G  {note}")
        L.append("-" * 96)
        L.append(f"{'TOTAL Acacia (object storage)':<{w1}}{'':>9}{'':>12}"
                 f"{arch_t/GB:>11.0f}G")
        L.append("")
        for it, arch, m, t, note in scratch:
            L.append(f"{it[:w1-1]:<{w1}}{'Scratch':>9}{m/GB:>11.1f}G"
                     f"{t/GB:>11.1f}G  {note}")
        L += ["=" * 96,
              "",
              "Acacia holds what must survive the project: reads, BAMs, gVCFs,",
              "call sets. Scratch holds Nextflow's work directory, which is",
              "transient and should be deleted once outputs are published --",
              "it is roughly the same size again while a run is live.",
              "",
              "The cohort VCF is scaled sub-linearly (square root): variant",
              "sites saturate as samples are added, while each added sample",
              "contributes one more genotype column per site. This is a",
              "conservative estimate, not a measurement -- confirm it against",
              "the full-cohort VCF once that exists.",
              "=" * 96]
        text = "\n".join(L)

    if args.output:
        open(args.output, "w").write(text + "\n")
        print(f"wrote {args.output}", file=sys.stderr)
    else:
        print(text)


if __name__ == "__main__":
    main()
