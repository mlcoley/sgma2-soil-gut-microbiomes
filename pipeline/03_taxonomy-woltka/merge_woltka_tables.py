#!/usr/bin/env python
"""Merge per-sample Woltka outputs into cohort matrices.

Per-sample inputs: <woltka>/<sample>/<rank>.biom (8 ranks) and
<woltka>/<sample>_processed.cov (genome_id, coverage).

Outputs in --out: cohort_<rank>.tsv/.biom, cohort_coverage.tsv (genome x sample),
cohort_classified_summary.tsv (per sample: compartment, post-host reads, classified
reads, % classified) and merge_report.txt. The sample list comes from the
read-accounting table; samples without output are reported as missing.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import biom
from biom.util import biom_open

RANKS = ["phylum", "class", "order", "family", "genus", "species", "none", "free"]


def load_biom_column(path, sample):
    """Load a single-sample .biom -> Series indexed by feature id, named `sample`."""
    t = biom.load_table(path)
    df = t.to_dataframe(dense=True)            # observations x samples
    if df.shape[1] != 1:
        # woltka --no-demux writes one sample; if not, sum columns to be safe
        s = df.sum(axis=1)
    else:
        s = df.iloc[:, 0]
    s.name = sample
    return s.astype(np.int64)


def write_biom(df, path, log):
    """df: features x samples (ints) -> BIOM 2.1 HDF5, fall back to JSON."""
    table = biom.Table(
        df.values,
        observation_ids=[str(x) for x in df.index],
        sample_ids=[str(x) for x in df.columns],
    )
    try:
        with biom_open(path, "w") as fh:
            table.to_hdf5(fh, "merge_woltka_tables.py")
    except Exception as e:  # h5py missing or other -> JSON fallback
        jpath = path + ".json"
        with open(jpath, "w") as fh:
            fh.write(table.to_json("merge_woltka_tables.py"))
        log.append(f"  [warn] HDF5 BIOM failed ({e}); wrote JSON BIOM -> {os.path.basename(jpath)}")


def write_tsv(df, path):
    out = df.copy()
    out.index.name = "#FeatureID"
    out.to_csv(path, sep="\t")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--woltka-dir", required=True)
    ap.add_argument("--read-accounting", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ranks", default=",".join(RANKS),
                    help="comma-separated subset of ranks to merge")
    args = ap.parse_args()

    ranks = [r.strip() for r in args.ranks.split(",") if r.strip()]
    os.makedirs(args.out, exist_ok=True)
    log = []

    # ---- authoritative sample list + metadata (compartment, posthost_reads) ----
    acc = pd.read_csv(args.read_accounting, sep="\t", dtype={"sample_id": str})
    acc = acc.set_index("sample_id")
    all_samples = list(acc.index)
    log.append(f"read-accounting samples: {len(all_samples)}")

    # samples that actually produced a per-sample dir (use species.biom as the sentinel)
    present = [s for s in all_samples
              if os.path.isfile(os.path.join(args.woltka_dir, s, "species.biom"))]
    missing = [s for s in all_samples if s not in present]
    log.append(f"samples with woltka output: {len(present)}")
    if missing:
        log.append(f"MISSING ({len(missing)}): {', '.join(missing)}")
        log.append("  (near-empty blanks/controls may legitimately produce no table; "
                   "check logs for any genuine failure)")

    # ---- merge each rank ----
    for rank in ranks:
        cols = {}
        skipped = []
        for s in present:
            p = os.path.join(args.woltka_dir, s, f"{rank}.biom")
            if not os.path.isfile(p):
                skipped.append(s)
                continue
            cols[s] = load_biom_column(p, s)
        if not cols:
            log.append(f"[{rank}] no tables found — skipped")
            continue
        df = pd.DataFrame(cols).reindex(columns=[s for s in present if s in cols])
        df = df.fillna(0).astype(np.int64)
        df = df.loc[df.sum(axis=1) > 0]  # drop all-zero features
        write_tsv(df, os.path.join(args.out, f"cohort_{rank}.tsv"))
        write_biom(df, os.path.join(args.out, f"cohort_{rank}.biom"), log)
        msg = f"[{rank}] {df.shape[0]} features x {df.shape[1]} samples"
        if skipped:
            msg += f"  (no {rank}.biom for {len(skipped)}: {', '.join(skipped)})"
        log.append(msg)

    # ---- merge per-genome coverage ----
    cov_cols = {}
    for p in sorted(glob.glob(os.path.join(args.woltka_dir, "*_processed.cov"))):
        s = os.path.basename(p)[:-len("_processed.cov")]
        ser = pd.read_csv(p, sep="\t", header=None, index_col=0).iloc[:, 0]
        ser.name = s
        cov_cols[s] = ser
    if cov_cols:
        order = [s for s in all_samples if s in cov_cols] + \
                [s for s in cov_cols if s not in set(all_samples)]
        cov = pd.DataFrame(cov_cols).reindex(columns=order).fillna(0.0)
        cov = cov.loc[cov.sum(axis=1) > 0]
        cov.index.name = "#GenomeID"
        cov.to_csv(os.path.join(args.out, "cohort_coverage.tsv"), sep="\t")
        log.append(f"[coverage] {cov.shape[0]} genomes x {cov.shape[1]} samples")
    else:
        log.append("[coverage] no *_processed.cov found")

    # ---- %classified per sample (classified = column sum of the per-genome 'none' rank) ----
    none_tsv = os.path.join(args.out, "cohort_none.tsv")
    if os.path.isfile(none_tsv):
        none = pd.read_csv(none_tsv, sep="\t", index_col=0)
        classified = none.sum(axis=0).astype(np.int64)
        summ = acc.loc[acc.index.isin(classified.index),
                       ["compartment", "posthost_reads"]].copy()
        summ["classified_reads"] = classified.reindex(summ.index).fillna(0).astype(np.int64)
        summ["pct_classified"] = (100.0 * summ["classified_reads"] /
                                  summ["posthost_reads"].replace(0, np.nan)).round(4)
        summ.to_csv(os.path.join(args.out, "cohort_classified_summary.tsv"), sep="\t")
        log.append("\n%classified by compartment (median, n):")
        for comp, g in summ.groupby("compartment"):
            med = g["pct_classified"].median(skipna=True)
            log.append(f"  {comp:>16}: median {med:6.2f}%   (n={len(g)})")
    else:
        log.append("[summary] cohort_none.tsv not produced — skipped %classified")

    report = "\n".join(log)
    with open(os.path.join(args.out, "merge_report.txt"), "w") as fh:
        fh.write(report + "\n")
    print(report)
    print(f"\nOutputs written to: {args.out}")


if __name__ == "__main__":
    sys.exit(main())
