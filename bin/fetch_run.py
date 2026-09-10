#!/usr/bin/env python3
"""
Download one FASTQ and verify its MD5.

Used by the SRA_FETCH process. Deliberately stdlib-only so it runs in any
Python container without needing curl or wget installed.

    fetch_run.py <url> <md5> <output>

Resumes a partial file with an HTTP Range request. A file that fails checksum
is deleted and the exit status is non-zero, so Nextflow's errorStrategy
retries it rather than passing a corrupt FASTQ to the aligner.
"""
import hashlib
import os
import sys
import time
import urllib.error
import urllib.request

CHUNK = 1 << 20  # 1 MiB


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def download(url, out, attempt):
    # Resume if we already have part of the file.
    start = os.path.getsize(out) if os.path.exists(out) else 0
    req = urllib.request.Request(url)
    mode = "wb"
    if start:
        req.add_header("Range", f"bytes={start}-")
        mode = "ab"

    try:
        with urllib.request.urlopen(req, timeout=300) as resp:
            # 200 to a Range request means the server ignored it: start over.
            if start and resp.status == 200:
                start, mode = 0, "wb"
            with open(out, mode) as fh:
                while True:
                    block = resp.read(CHUNK)
                    if not block:
                        break
                    fh.write(block)
    except urllib.error.HTTPError as exc:
        # 416 means we already have the whole file.
        if exc.code == 416 and start:
            return
        raise


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    url, want, out = sys.argv[1:4]

    for attempt in range(1, 6):
        try:
            if not (os.path.exists(out) and md5_of(out) == want):
                download(url, out, attempt)
        except Exception as exc:                              # noqa: BLE001
            print(f"attempt {attempt}: transfer failed: {exc}", file=sys.stderr)
            time.sleep(min(60, 10 * attempt))
            continue

        if not os.path.exists(out):
            print(f"attempt {attempt}: no output produced", file=sys.stderr)
            continue

        got = md5_of(out)
        if got == want:
            print(f"OK {out} ({os.path.getsize(out):,} bytes)")
            return 0

        print(f"attempt {attempt}: MD5 mismatch for {out} "
              f"(expected {want}, got {got}) -- refetching from scratch",
              file=sys.stderr)
        os.remove(out)
        time.sleep(min(60, 10 * attempt))

    sys.exit(f"ERROR: failed to download {url} with matching MD5 after 5 attempts")


if __name__ == "__main__":
    sys.exit(main())
