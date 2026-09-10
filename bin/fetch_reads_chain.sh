#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# fetch_reads_chain.sh -- download the read set as a chain of short SLURM jobs.
#
# A single long download job is effectively unschedulable on a busy cluster:
# it is too large to backfill, so it waits for a full-priority slot. Because
# fetch_reads.sh is resumable and checksum-verified, the download can instead
# be split into many short, single-core jobs chained with --dependency=afterany.
# Each is small enough to slot into a backfill gap, and each resumes exactly
# where the previous one stopped.
#
#   bin/fetch_reads_chain.sh /group/peg/cicer/cret/reads [n_jobs] [hours] [par]
#
# The chain stops early on its own: once every file verifies, each remaining
# job exits within seconds.
# ---------------------------------------------------------------------------
set -euo pipefail

OUTDIR=${1:?usage: fetch_reads_chain.sh <outdir> [n_jobs] [hours] [parallel]}
NJOBS=${2:-12}
HOURS=${3:-3}
PAR=${4:-8}
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$PROJECT_DIR/logs"

dep=""
for i in $(seq 1 "$NJOBS"); do
    jid=$(sbatch --parsable \
        ${dep:+--dependency=afterany:$dep} \
        --job-name="Crfetch$i" \
        --partition=work \
        --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=2G \
        --time="${HOURS}:00:00" \
        --output="$PROJECT_DIR/logs/Cr_fetch_chain_%j.out" \
        --error="$PROJECT_DIR/logs/Cr_fetch_chain_%j.err" \
        --wrap "cd '$PROJECT_DIR' && bin/fetch_reads.sh -o '$OUTDIR' -j $PAR")
    echo "  link $i/$NJOBS -> job $jid${dep:+ (after $dep)}"
    dep=$jid
done
echo
echo "Chain submitted. Monitor with:  squeue -u \$USER"
echo "Verify at any time with:        bin/fetch_reads.sh -o '$OUTDIR' -c"
