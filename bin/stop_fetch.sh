#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# stop_fetch.sh -- stop a running fetch_reads.sh cleanly.
#
# fetch_reads.sh runs under `setsid`, so the whole job is one process group and
# can be signalled atomically. Killing only the parent leaves the xargs layer
# and its curl children running, which is what allows two fetchers to overlap.
#
#   bin/stop_fetch.sh /group/peg/cicer/cret/reads
# ---------------------------------------------------------------------------
set -uo pipefail
OUTDIR=${1:?usage: stop_fetch.sh <outdir>}
PGIDFILE="$OUTDIR/.fetch.pgid"

if [[ -f "$PGIDFILE" ]]; then
    pgid=$(cat "$PGIDFILE")
    if kill -0 "-$pgid" 2>/dev/null; then
        echo "stopping process group $pgid ..."
        kill -TERM "-$pgid" 2>/dev/null; sleep 5
        kill -KILL "-$pgid" 2>/dev/null; sleep 2
    fi
    rm -f "$PGIDFILE"
fi

# Belt and braces: any curl still pointed at ENA, regardless of group.
for round in 1 2 3; do
    n=0
    for pid in $(pgrep -x curl 2>/dev/null); do
        if tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null | grep -q 'sra.ebi.ac.uk'; then
            kill -9 "$pid" 2>/dev/null && n=$((n+1))
        fi
    done
    [[ $n -eq 0 ]] && break
    echo "  killed $n stray curl"
    sleep 2
done

rm -f "$OUTDIR/.fetch.lock"
echo "stopped. Re-run bin/fetch_reads.sh to resume; partial files are resumed"
echo "and any file failing MD5 is refetched from scratch."
