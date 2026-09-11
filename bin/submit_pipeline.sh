#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# submit_pipeline.sh -- submit the Nextflow driver with an escalating-memory
# retry ladder.
#
# Two independent things can run out of memory:
#
#   individual tasks   handled inside the pipeline. conf/base.config scales
#                      memory with task.attempt and retries the OOM/kill exit
#                      codes (104,134,137,139,140,143,247) up to three times.
#
#   the driver itself  handled here. The Nextflow driver is a long-lived JVM
#                      whose footprint grows with the number of tasks it is
#                      tracking, so a cohort of 161 samples x 8 intervals needs
#                      more than a 2-sample test. A driver OOM kills the whole
#                      run regardless of how well-sized the tasks were.
#
# This submits a chain of attempts with doubling memory, each depending on the
# previous with --dependency=afternotok, so an attempt runs only if the one
# before it failed. Every attempt passes -resume, so a retry continues from
# cached work rather than starting over. If an attempt succeeds, the remaining
# links never run.
#
#   bin/submit_pipeline.sh -- -profile uwa,apptainer --fasta ... --outdir ...
#
#   -m <MB>      memory for the first attempt   (default 8192)
#   -n <int>     number of attempts             (default 4)
#   -t <time>    wall clock per attempt         (default 24:00:00)
#   -c <int>     cpus for the driver            (default 2)
#   --           everything after is passed to `nextflow run`
# ---------------------------------------------------------------------------
set -euo pipefail

MEM=8192; ATTEMPTS=4; WALL="24:00:00"; CPUS=2
while [[ $# -gt 0 ]]; do
    case $1 in
        -m) MEM=$2; shift 2 ;;
        -n) ATTEMPTS=$2; shift 2 ;;
        -t) WALL=$2; shift 2 ;;
        -c) CPUS=$2; shift 2 ;;
        --) shift; break ;;
        -h|--help) sed -n '2,30p' "$0"; exit 0 ;;
        *)  break ;;
    esac
done

[[ $# -eq 0 ]] && { echo "ERROR: no nextflow arguments given (use -- before them)" >&2; exit 2; }

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$PROJECT_DIR/logs"

# -resume is what makes escalation cheap: each retry reuses everything the
# previous attempt completed. Added once here if the caller has not.
NF_ARGS=("$@")
[[ " ${NF_ARGS[*]} " == *" -resume "* ]] || NF_ARGS+=(-resume)

echo "driver retry ladder: $ATTEMPTS attempts, ${MEM}MB doubling, ${WALL} each"
echo

dep=""
mem=$MEM
for i in $(seq 1 "$ATTEMPTS"); do
    jid=$(sbatch --parsable \
        ${dep:+--dependency=afternotok:$dep} \
        --job-name="PEalign${i}" \
        --partition=work \
        --nodes=1 --ntasks=1 --cpus-per-task="$CPUS" \
        --mem="${mem}M" --time="$WALL" \
        --output="$PROJECT_DIR/logs/PE_align_a${i}_%j.out" \
        --error="$PROJECT_DIR/logs/PE_align_a${i}_%j.err" \
        "$PROJECT_DIR/bin/run_pipeline.sbatch" "${NF_ARGS[@]}")
    if [[ -n "$dep" ]]; then
        printf "  attempt %d: job %-9s mem=%6sMB   (runs only if %s fails)\n" \
               "$i" "$jid" "$mem" "$dep"
    else
        printf "  attempt %d: job %-9s mem=%6sMB   (runs immediately)\n" \
               "$i" "$jid" "$mem"
    fi
    dep=$jid
    mem=$(( mem * 2 ))
done

echo
echo "Each attempt runs only if the previous one FAILED (--dependency=afternotok)."
echo "Cancel the whole ladder with:  scancel --name=PEalign1 --name=PEalign2 ..."
echo "or:  squeue -u \$USER -h -o '%i %j' | awk '\$2 ~ /^PEalign/ {print \$1}' | xargs -r scancel"
