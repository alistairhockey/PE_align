#!/bin/bash
#SBATCH --job-name=PE_align
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/PE_align_%j.out
#SBATCH --error=logs/PE_align_%j.err
# ---------------------------------------------------------------------------
# Submit the Nextflow head job.
#
# Nothing in this pipeline should run on a login node. This script runs the
# Nextflow *driver* inside a small allocation; the driver then submits each
# task as its own SLURM job via the `slurm` executor in the chosen profile.
# The driver itself is light -- 2 cores and 8 GB is ample -- but it must live
# for the whole run, hence the long wall time.
#
#   sbatch bin/run_pipeline.sbatch -profile uwa,apptainer \
#       --sra_metadata assets/cret_metadata.tsv \
#       --fasta /path/to/assembly.fna \
#       --reads_dir /group/peg/cicer/cret/reads \
#       --outdir results --benchmark_label full_v1 -resume
#
# Everything after the script name is passed straight to `nextflow run`.
# ---------------------------------------------------------------------------
set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$PROJECT_DIR"
mkdir -p logs

# Keep container and Nextflow caches off /home, which is small and often
# not visible to compute nodes in the same way /group is.
export NXF_APPTAINER_CACHEDIR="${NXF_APPTAINER_CACHEDIR:-$PROJECT_DIR/.apptainer_cache}"
export NXF_SINGULARITY_CACHEDIR="$NXF_APPTAINER_CACHEDIR"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-${TMPDIR:-/tmp}}"
export NXF_OPTS="${NXF_OPTS:--Xms1g -Xmx6g}"
mkdir -p "$NXF_APPTAINER_CACHEDIR"

echo "host        : $(hostname)"
echo "job         : ${SLURM_JOB_ID:-none}"
echo "project     : $PROJECT_DIR"
echo "cache       : $NXF_APPTAINER_CACHEDIR"
echo "nextflow    : $(nextflow -version 2>&1 | sed -n 's/.*version \(.*\)/\1/p' | head -1)"
echo "args        : $*"
echo

nextflow run "$PROJECT_DIR" "$@"
