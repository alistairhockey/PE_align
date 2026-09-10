#!/usr/bin/env python3
"""
Turn an NCBI SRA run table into a normalised run manifest.

Accepts whatever the SRA Run Selector hands you: UTF-8 or UTF-16, comma or tab
separated, with or without a BOM, and with the column-name variations NCBI uses
across its export paths. Resolves each run against the ENA file report to get
paired FASTQ URLs, MD5 checksums and byte counts.

    parse_runtable.py SraRunTable.txt -o runs.tsv [--layout PAIRED] [--assay WGS]

Output columns:
    run  sample  population  biosample  bases  bytes
    url_1  md5_1  bytes_1  url_2  md5_2  bytes_2
"""
import argparse
import csv
import io
import json
import sys
import urllib.parse
import urllib.request

ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = ("run_accession,fastq_ftp,fastq_md5,fastq_bytes,"
              "library_layout,scientific_name,sample_alias,base_count")

# NCBI is not consistent about these; accept any of them.
RUN_KEYS    = ("Run", "run", "run_accession", "Run accession")
SAMPLE_KEYS = ("Sample Name", "SampleName", "sample_name", "isolate",
               "Library Name", "sample_alias", "BioSample")
ASSAY_KEYS  = ("Assay Type", "assay_type", "LibraryStrategy")
LAYOUT_KEYS = ("LibraryLayout", "Library Layout", "library_layout")
BASES_KEYS  = ("Bases", "bases", "base_count")
BYTES_KEYS  = ("Bytes", "bytes")
BIOS_KEYS   = ("BioSample", "biosample", "BioSample Accession")
STUDY_KEYS  = ("BioProject", "SRA Study", "study_accession", "bioproject",
               "secondary_study_accession")


def read_table(path):
    """Decode the file regardless of encoding, then sniff the delimiter."""
    raw = open(path, "rb").read()
    for enc in ("utf-16", "utf-8-sig", "utf-8", "latin-1"):
        try:
            text = raw.decode(enc)
            # A mis-guessed encoding usually shows up as interleaved NULs.
            if "\x00" in text:
                continue
            break
        except (UnicodeDecodeError, UnicodeError):
            continue
    else:
        sys.exit(f"ERROR: could not decode {path} as UTF-16, UTF-8 or Latin-1")

    text = text.replace("\r\n", "\n").replace("\r", "\n")
    header = text.split("\n", 1)[0]
    delim = "\t" if header.count("\t") >= header.count(",") else ","
    rows = list(csv.DictReader(io.StringIO(text), delimiter=delim))
    if not rows:
        sys.exit(f"ERROR: no data rows parsed from {path}")
    return rows


def pick(row, keys, default=""):
    for k in keys:
        v = row.get(k)
        if v not in (None, ""):
            return v.strip()
    return default


def _ena_get(params, timeout=180, attempts=4):
    """GET an ENA portal endpoint, returning parsed TSV rows."""
    endpoint = params.pop("_endpoint")
    url = f"https://www.ebi.ac.uk/ena/portal/api/{endpoint}?" + urllib.parse.urlencode(params)
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as fh:
                text = fh.read().decode("utf-8")
            return list(csv.DictReader(io.StringIO(text), delimiter="\t"))
        except Exception as exc:                              # noqa: BLE001
            if attempt == attempts:
                sys.exit(f"ERROR: ENA request failed after {attempts} attempts: {exc}")
            print(f"  ENA retry {attempt} ({exc})", file=sys.stderr)
    return []


def ena_lookup(accessions, studies=(), chunk=100):
    """
    Resolve run accessions to FASTQ URLs, MD5s and sizes.

    Two strategies, because the portal API's `accession=` parameter silently
    returns an empty result set for a comma-separated list -- it accepts one
    accession, or one study/project:

      1. If the run table names a study or BioProject, fetch that study once
         and keep the runs we care about. One request for any cohort size.
      2. Otherwise fall back to the `search` endpoint with batched
         `run_accession="X" OR ...` queries.
    """
    want = set(accessions)
    out = {}

    for study in studies:
        rows = _ena_get({"_endpoint": "filereport", "accession": study,
                         "result": "read_run", "fields": ENA_FIELDS,
                         "format": "tsv", "limit": "0"})
        for r in rows:
            if r.get("run_accession") in want:
                out[r["run_accession"]] = r
        print(f"  study {study}: {len(rows)} runs, {len(out)}/{len(want)} matched",
              file=sys.stderr)
        if len(out) == len(want):
            return out

    missing = sorted(want - set(out))
    if missing:
        print(f"  resolving {len(missing)} run(s) individually", file=sys.stderr)
    for i in range(0, len(missing), chunk):
        batch = missing[i:i + chunk]
        q = " OR ".join(f'run_accession="{a}"' for a in batch)
        rows = _ena_get({"_endpoint": "search", "result": "read_run",
                         "query": q, "fields": ENA_FIELDS,
                         "format": "tsv", "limit": "0"})
        for r in rows:
            out[r["run_accession"]] = r
        print(f"  resolved {len(out)}/{len(want)}", file=sys.stderr)

    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runtable")
    ap.add_argument("-o", "--output", default="runs.tsv")
    ap.add_argument("--layout", default="PAIRED",
                    help="keep only runs with this LibraryLayout ('ANY' to disable)")
    ap.add_argument("--assay", default=None,
                    help="keep only runs with this Assay Type, e.g. WGS")
    ap.add_argument("--organism", default=None,
                    help="keep only runs whose ENA scientific_name matches")
    ap.add_argument("--summary", default=None, help="write a JSON summary here")
    args = ap.parse_args()

    rows = read_table(args.runtable)
    print(f"read {len(rows)} rows from {args.runtable}", file=sys.stderr)

    if not any(k in rows[0] for k in RUN_KEYS):
        sys.exit(f"ERROR: no run-accession column found. Columns: "
                 f"{', '.join(list(rows[0])[:15])}")

    kept, dropped = [], {"layout": 0, "assay": 0, "no_run": 0}
    for r in rows:
        run = pick(r, RUN_KEYS)
        if not run:
            dropped["no_run"] += 1
            continue
        if args.assay:
            a = pick(r, ASSAY_KEYS)
            if a and a.upper() != args.assay.upper():
                dropped["assay"] += 1
                continue
        if args.layout.upper() != "ANY":
            lay = pick(r, LAYOUT_KEYS)
            if lay and lay.upper() != args.layout.upper():
                dropped["layout"] += 1
                continue
        kept.append(r)

    print(f"kept {len(kept)} runs after run-table filters "
          f"(dropped: {dropped})", file=sys.stderr)
    if not kept:
        sys.exit("ERROR: no runs survived filtering")

    print("querying ENA file report ...", file=sys.stderr)
    studies = sorted({pick(r, STUDY_KEYS) for r in kept} - {""})
    if studies:
        print(f"  study accessions in run table: {', '.join(studies)}", file=sys.stderr)
    ena = ena_lookup({pick(r, RUN_KEYS) for r in kept}, studies=studies)

    cols = ["run", "sample", "population", "biosample", "bases", "bytes",
            "url_1", "md5_1", "bytes_1", "url_2", "md5_2", "bytes_2"]
    written, skipped = 0, {"not_in_ena": [], "not_paired": [], "organism": []}

    with open(args.output, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in sorted(kept, key=lambda x: pick(x, RUN_KEYS)):
            run = pick(r, RUN_KEYS)
            e = ena.get(run)
            if not e:
                skipped["not_in_ena"].append(run)
                continue
            if args.organism and args.organism.lower() not in e.get(
                    "scientific_name", "").lower():
                skipped["organism"].append(run)
                continue

            urls = [u for u in e.get("fastq_ftp", "").split(";") if u]
            md5s = [m for m in e.get("fastq_md5", "").split(";") if m]
            byts = [b for b in e.get("fastq_bytes", "").split(";") if b]

            # ENA sometimes lists an unpaired "orphan" file alongside R1/R2;
            # select the true mates by suffix rather than by position.
            pair = {}
            for u, m, b in zip(urls, md5s, byts):
                if u.endswith("_1.fastq.gz"):
                    pair[1] = (u, m, b)
                elif u.endswith("_2.fastq.gz"):
                    pair[2] = (u, m, b)
            if 1 not in pair or 2 not in pair:
                skipped["not_paired"].append(run)
                continue

            sample = pick(r, SAMPLE_KEYS) or e.get("sample_alias") or run
            sample = sample.replace(" ", "_")
            w.writerow([
                run, sample, sample.rsplit("_", 1)[0],
                pick(r, BIOS_KEYS),
                pick(r, BASES_KEYS) or e.get("base_count", ""),
                pick(r, BYTES_KEYS),
                "https://" + pair[1][0], pair[1][1], pair[1][2],
                "https://" + pair[2][0], pair[2][1], pair[2][2],
            ])
            written += 1

    for reason, runs in skipped.items():
        if runs:
            print(f"WARNING: skipped {len(runs)} run(s) [{reason}]: "
                  f"{', '.join(runs[:5])}{' ...' if len(runs) > 5 else ''}",
                  file=sys.stderr)

    print(f"wrote {written} runs to {args.output}", file=sys.stderr)
    if written == 0:
        sys.exit("ERROR: no runs resolved to paired FASTQs")

    if args.summary:
        json.dump({"runs": written,
                   "skipped": {k: len(v) for k, v in skipped.items()},
                   "dropped": dropped},
                  open(args.summary, "w"), indent=2)


if __name__ == "__main__":
    main()
