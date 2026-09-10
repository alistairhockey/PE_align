#!/usr/bin/env python3
"""
Generate a self-contained test dataset.

Creates a small synthetic reference and simulates paired-end reads from it for
several samples, one of which has two sequencing runs so the per-sample merge
path is exercised. No external data required, so the test runs anywhere.

    tests/make_test_data.py -o tests/data
"""
import argparse
import gzip
import os
import random

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def make_reference(path, contigs, length, seed):
    rng = random.Random(seed)
    seqs = {}
    with open(path, "w") as fh:
        for i in range(1, contigs + 1):
            name = f"testchr{i}"
            # Skew composition slightly so the sequence is not uniformly random,
            # which keeps bwa's seeding behaviour closer to a real genome.
            seq = "".join(rng.choices("ACGT", weights=[3, 2, 2, 3], k=length))
            seqs[name] = seq
            fh.write(f">{name}\n")
            for j in range(0, length, 60):
                fh.write(seq[j:j + 60] + "\n")
    return seqs


def simulate(seqs, out1, out2, n_pairs, read_len, frag, err, snp_rate, seed):
    rng = random.Random(seed)
    names = list(seqs)

    # Give each sample its own variants so the cohort is polymorphic and the
    # joint-genotyping and PCA steps have something to work with.
    variants = {}
    for name, seq in seqs.items():
        variants[name] = {
            pos: rng.choice([b for b in "ACGT" if b != seq[pos]])
            for pos in rng.sample(range(len(seq)), int(len(seq) * snp_rate))
        }

    with gzip.open(out1, "wt") as f1, gzip.open(out2, "wt") as f2:
        for i in range(n_pairs):
            name = rng.choice(names)
            seq = seqs[name]
            start = rng.randrange(0, len(seq) - frag)
            frag_seq = list(seq[start:start + frag])

            for pos, alt in variants[name].items():
                if start <= pos < start + frag:
                    frag_seq[pos - start] = alt

            r1 = "".join(frag_seq[:read_len])
            r2 = rc("".join(frag_seq[-read_len:]))

            def sequencing_error(read):
                out = []
                for base in read:
                    if rng.random() < err:
                        out.append(rng.choice([b for b in "ACGT" if b != base]))
                    else:
                        out.append(base)
                return "".join(out)

            r1, r2 = sequencing_error(r1), sequencing_error(r2)
            q = "I" * read_len
            f1.write(f"@read{i}/1\n{r1}\n+\n{q}\n")
            f2.write(f"@read{i}/2\n{r2}\n+\n{q}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", default="tests/data")
    ap.add_argument("--contigs", type=int, default=2)
    ap.add_argument("--length", type=int, default=120_000)
    ap.add_argument("--pairs", type=int, default=12_000)
    ap.add_argument("--read-len", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ref = os.path.join(args.outdir, "test_ref.fna")
    seqs = make_reference(ref, args.contigs, args.length, args.seed)
    print(f"reference: {ref} ({args.contigs} x {args.length:,} bp)")

    # SAMPLE_A has two runs -- this is what exercises SAMTOOLS_MERGE.
    plan = [("SAMPLE_A", "run1"), ("SAMPLE_A", "run2"),
            ("SAMPLE_B", "run1"), ("SAMPLE_C", "run1")]

    rows = ["sample,run,fastq_1,fastq_2"]
    for idx, (sample, run) in enumerate(plan):
        o1 = os.path.abspath(os.path.join(args.outdir, f"{run}_{sample}_1.fastq.gz"))
        o2 = os.path.abspath(os.path.join(args.outdir, f"{run}_{sample}_2.fastq.gz"))
        simulate(seqs, o1, o2, args.pairs, args.read_len,
                 frag=args.read_len * 3, err=0.002, snp_rate=0.001,
                 seed=args.seed + 1000 * (ord(sample[-1]) - 65) + idx)
        rows.append(f"{sample},{run},{o1},{o2}")
        print(f"  {sample}/{run}: {args.pairs:,} pairs")

    sheet = os.path.join(args.outdir, "..", "test_samplesheet.csv")
    sheet = os.path.normpath(sheet)
    open(sheet, "w").write("\n".join(rows) + "\n")
    print(f"samplesheet: {sheet}")


if __name__ == "__main__":
    main()
