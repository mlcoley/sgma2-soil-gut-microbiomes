#!/usr/bin/env python3
"""Pull abundance rows for the selected background genes.

Inputs: --abundance gene_abundance_long.tsv.gz (no header), --genes (gene ids, one per
line), --cohort T23c_cohort_compartments.tsv (analysed samples).
Output: --out TSV with gene_id, sample_id, count, covered_fraction (cohort samples
only). A header row in the abundance table is detected, not assumed.
"""

import argparse, gzip, os
from itertools import chain


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def log(msg):
    print(msg, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--abundance", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-cohort", type=int, default=194)
    a = ap.parse_args()

    log("== background abundance extraction ==")

    with open(a.cohort) as fh:
        fh.readline()
        cohort = {ln.split("\t")[0] for ln in fh if ln.strip()}
    log(f"cohort samples: {len(cohort)}")
    assert len(cohort) == a.n_cohort, f"expected {a.n_cohort} cohort samples, got {len(cohort)}"

    with open(a.genes) as fh:
        want = {ln.strip() for ln in fh if ln.strip()}
    log(f"background genes wanted: {len(want):,}")
    assert want, "empty gene list"

    tmp = a.out + ".partial"
    n_rows = n_kept = 0
    seen = set()
    with opener(a.abundance) as fh, gzip.open(tmp, "wt") as out:
        out.write("gene_id\tsample_id\tcount\tcovered_fraction\n")

        first = fh.readline()
        p0 = first.rstrip("\n").split("\t")
        if p0[:1] == ["gene_id"]:
            log("  header row detected and consumed")
            stream = fh
        else:
            assert len(p0) >= 4 and p0[2].isdigit(), \
                f"first line is neither a header nor a 4-column data row: {p0}"
            log("  no header row — treating line 1 as data (expected for this file)")
            stream = chain([first], fh)

        for ln in stream:
            n_rows += 1
            if n_rows % 50_000_000 == 0:
                log(f"  rows {n_rows:,}  kept {n_kept:,}  genes seen {len(seen):,}")
            p = ln.split("\t", 2)
            if p[0] not in want or p[1] not in cohort:
                continue
            out.write(ln if ln.endswith("\n") else ln + "\n")
            seen.add(p[0])
            n_kept += 1

    os.replace(tmp, a.out)
    log(f"rows read {n_rows:,}, kept {n_kept:,}")
    log(f"genes with >=1 cohort row: {len(seen):,} of {len(want):,} wanted")

    missing = len(want) - len(seen)
    if missing:
        # Selected genes were detected in the cohort, so a gene with no rows here
        # means the inputs disagree.
        raise SystemExit(
            f"ERROR: {missing:,} selected genes returned no cohort rows; the gene "
            f"list and the abundance table are out of sync.")
    assert n_kept > 0, "no rows written — refusing to publish an empty extraction"
    log("done.")


if __name__ == "__main__":
    main()
