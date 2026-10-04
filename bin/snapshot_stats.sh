#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# snapshot_stats.sh -- freeze this run's resource statistics before they are lost.
#
# SLURM accounting is rotated, and this cluster's job-ID counter has already
# reset once, which truncates what sacct will return. Nextflow traces can also
# be lost (a pre-existing trace file disables the writer). A long project
# therefore needs its statistics captured as they happen, not reconstructed
# afterwards.
#
# Writes, into stats/<label>/:
#   sacct-raw.psv     every accounting record, pipe-separated and unprocessed
#   summary.txt       per-process CPU and memory summary
#   trace-*.txt       any Nextflow traces found under --results
#   MANIFEST.txt      what was captured, when, and from where
#
#   bin/snapshot_stats.sh -l alignment -s 2026-09-10
#   bin/snapshot_stats.sh -l genotyping -s 2026-10-04 -r results/from_bams
# ---------------------------------------------------------------------------
set -uo pipefail

LABEL=""; SINCE=""; RESULTS=""; OUTROOT="stats"
while getopts "l:s:r:o:h" o; do case $o in
  l) LABEL=$OPTARG ;; s) SINCE=$OPTARG ;; r) RESULTS=$OPTARG ;;
  o) OUTROOT=$OPTARG ;; h) sed -n '2,20p' "$0"; exit 0 ;;
  *) exit 2 ;;
esac; done

[[ -z "$LABEL" ]] && { echo "ERROR: -l <label> required" >&2; exit 2; }
[[ -z "$SINCE" ]] && { echo "ERROR: -s <since date> required" >&2; exit 2; }

HERE="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$HERE/$OUTROOT/$LABEL"
mkdir -p "$OUT"

FIELDS="JobID,JobName%200,State,Elapsed,TotalCPU,AllocCPUS,ReqMem,MaxRSS,MaxVMSize,MaxDiskRead,MaxDiskWrite,Start"

echo "capturing accounting since $SINCE ..."
sacct -u "$USER" -S "$SINCE" -P -n --format="$FIELDS" > "$OUT/sacct-raw.psv" 2>/dev/null
n=$(wc -l < "$OUT/sacct-raw.psv")
nf=$(grep -c '|nf-' "$OUT/sacct-raw.psv" || true)
echo "  $n records, $nf nextflow task records"

if [[ "${nf:-0}" -eq 0 ]]; then
    echo "WARNING: no nf- records captured. Either the window is wrong, or" >&2
    echo "         accounting has already rotated past this run." >&2
fi

echo "summarising ..."
"$HERE/bin/summarise_sacct.py" --from-file "$OUT/sacct-raw.psv" \
    > "$OUT/summary.txt" 2>/dev/null || echo "  (summary failed)"

if [[ -n "$RESULTS" ]]; then
    find "$RESULTS" -name "trace-*.txt" -size +1k 2>/dev/null | while read -r t; do
        cp "$t" "$OUT/$(basename "$t")" && echo "  kept trace $(basename "$t")"
    done
fi

{
  echo "label:        $LABEL"
  echo "captured:     $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "sacct since:  $SINCE"
  echo "host:         $(hostname)"
  echo "user:         $USER"
  echo "git commit:   $(cd "$HERE" && git rev-parse --short HEAD 2>/dev/null)"
  echo "records:      $n total, $nf nextflow tasks"
  echo "traces kept:  $(ls "$OUT"/trace-*.txt 2>/dev/null | wc -l)"
  echo
  echo "Combine snapshots later with:"
  echo "  bin/allocation_table.py $(for d in "$HERE/$OUTROOT"/*/; do printf -- '--from-file %s/sacct-raw.psv ' "$d"; done)"
} > "$OUT/MANIFEST.txt"

echo
cat "$OUT/MANIFEST.txt"
echo
echo "snapshot written to $OUT"
