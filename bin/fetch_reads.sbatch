#!/bin/bash
#SBATCH --job-name=Cr_fetch
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/Cr_fetch_%j.out
#SBATCH --error=logs/Cr_fetch_%j.err
# ---------------------------------------------------------------------------
# Bulk-download the read set on a compute node.
#
#   sbatch bin/fetch_reads.sbatch /group/peg/cicer/cret/reads [jobs]
#
# The work is network- and checksum-bound rather than CPU-bound, so a small
# allocation is right. Resumable: re-submitting after a timeout picks up where
# it left off, and every file is MD5-verified before being accepted.
# ---------------------------------------------------------------------------
set -euo pipefail

OUTDIR=${1:?usage: sbatch bin/fetch_reads.sbatch <outdir> [jobs]}
JOBS=${2:-8}
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$PROJECT_DIR"
mkdir -p logs

echo "host   : $(hostname)"
echo "outdir : $OUTDIR"
echo "jobs   : $JOBS"
echo

exec "$PROJECT_DIR/bin/fetch_reads.sh" -o "$OUTDIR" -j "$JOBS"
