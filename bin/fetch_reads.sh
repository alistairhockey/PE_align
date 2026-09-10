#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# fetch_reads.sh -- download the Cr_align read set from ENA.
#
# Resumable and idempotent: a file whose MD5 already matches the manifest is
# skipped, a partial file is resumed with `curl -C -`, and a file that fails
# checksum is deleted and retried. Safe to re-run after an interruption.
#
#   bin/fetch_reads.sh -o /group/peg/cicer/cret/reads [-j 6] [-m assets/cret_runs.tsv]
#
#   -o  output directory (required)
#   -j  parallel downloads (default 6; ENA throttles aggressively above ~8)
#   -m  manifest TSV (default assets/cret_runs.tsv)
#   -n  dry run: report what would be fetched, download nothing
#   -c  check only: verify existing files against the manifest and exit
# ---------------------------------------------------------------------------
set -uo pipefail

MANIFEST="$(dirname "$0")/../assets/cret_runs.tsv"
OUTDIR=""; JOBS=6; DRY=0; CHECKONLY=0

while getopts "o:j:m:nch" opt; do
  case $opt in
    o) OUTDIR=$OPTARG ;;
    j) JOBS=$OPTARG ;;
    m) MANIFEST=$OPTARG ;;
    n) DRY=1 ;;
    c) CHECKONLY=1 ;;
    h) sed -n '2,16p' "$0"; exit 0 ;;
    *) exit 2 ;;
  esac
done

[[ -z "$OUTDIR" ]] && { echo "ERROR: -o <outdir> is required" >&2; exit 2; }
[[ -f "$MANIFEST" ]] || { echo "ERROR: manifest not found: $MANIFEST" >&2; exit 2; }
mkdir -p "$OUTDIR" || exit 1

LOGDIR="$OUTDIR/.fetch_logs"; mkdir -p "$LOGDIR"

# One download unit: url, expected md5, destination filename.
fetch_one() {
  local url=$1 md5=$2 dest=$3 log="$LOGDIR/$(basename "$3").log"

  if [[ -f "$dest" ]]; then
    local have
    have=$(md5sum "$dest" | cut -d' ' -f1)
    if [[ "$have" == "$md5" ]]; then
      echo "OK-CACHED  $(basename "$dest")"
      return 0
    fi
    if [[ $CHECKONLY -eq 1 ]]; then
      echo "BAD-MD5    $(basename "$dest")  (expected $md5, got $have)"
      return 1
    fi
    # Size below expectation implies a truncated transfer -> resume.
    echo "RESUMING   $(basename "$dest")"
  fi

  [[ $CHECKONLY -eq 1 ]] && { echo "MISSING    $(basename "$dest")"; return 1; }
  [[ $DRY -eq 1 ]]       && { echo "WOULD-GET  $(basename "$dest")"; return 0; }

  local attempt
  for attempt in 1 2 3 4 5; do
    if curl -sSL --fail --retry 3 --retry-delay 10 --connect-timeout 30 \
            --speed-limit 10240 --speed-time 120 \
            -C - -o "$dest" "$url" >>"$log" 2>&1; then
      local have
      have=$(md5sum "$dest" | cut -d' ' -f1)
      if [[ "$have" == "$md5" ]]; then
        echo "OK         $(basename "$dest")"
        return 0
      fi
      echo "MD5-FAIL   $(basename "$dest") attempt $attempt -- refetching from scratch" >&2
      rm -f "$dest"
    else
      echo "CURL-FAIL  $(basename "$dest") attempt $attempt" >&2
      sleep $(( attempt * 15 ))
    fi
  done
  echo "FAILED     $(basename "$dest")" >&2
  return 1
}
export -f fetch_one
export LOGDIR DRY CHECKONLY

# Manifest -> one "url md5 dest" triple per line, R1 and R2 separately.
awk -F'\t' -v OUT="$OUTDIR" 'NR>1 {
    printf "%s\t%s\t%s/%s_%s_1.fastq.gz\n", $7,  $8, OUT, $1, $2;
    printf "%s\t%s\t%s/%s_%s_2.fastq.gz\n", $10, $11, OUT, $1, $2;
}' "$MANIFEST" \
| xargs -P "$JOBS" -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r u m d <<< "{}"; fetch_one "$u" "$m" "$d"'

echo
echo "---- summary ----"
echo "expected files : $(( ($(wc -l < "$MANIFEST") - 1) * 2 ))"
echo "present files  : $(find "$OUTDIR" -maxdepth 1 -name '*.fastq.gz' | wc -l)"
echo "total size     : $(du -sh "$OUTDIR" 2>/dev/null | cut -f1)"
