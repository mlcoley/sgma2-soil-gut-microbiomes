#!/usr/bin/env python
"""Per-sample MetaPhlAn classified summary and per-compartment medians.

Inputs: merged MetaPhlAn table (includes the UNCLASSIFIED row from
--unclassified_estimation) and the read-accounting TSV (compartment per sample).
Output: TSV with sample_id, compartment, pct_unclassified, pct_classified,
n_species. Prints median pct_classified per compartment.
Usage: metaphlan_classified_summary.py MERGED ACCOUNTING OUT
"""
import csv, sys, statistics as st

def sid(col):
    """Merged-table columns are '<P###>_metaphlan'; accounting keys on bare 'P###'."""
    return col[:-len("_metaphlan")] if col.endswith("_metaphlan") else col

def main():
    merged, acc_path, out = sys.argv[1], sys.argv[2], sys.argv[3]
    comp = {}
    with open(acc_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            comp[row["sample_id"]] = row["compartment"]

    rows = [r for r in csv.reader(open(merged), delimiter="\t") if r and not r[0].startswith("#mpa")]
    hdr = rows[0]
    samples = hdr[1:]
    unclass = {s: None for s in samples}
    nsp = {s: 0 for s in samples}
    for r in rows[1:]:
        clade = r[0]
        if clade == "UNCLASSIFIED":
            for i, s in enumerate(samples):
                unclass[s] = float(r[i + 1])
        elif clade.split("|")[-1].startswith("s__") and "t__" not in clade:
            for i, s in enumerate(samples):
                if float(r[i + 1]) > 0:
                    nsp[s] += 1

    with open(out, "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "compartment", "pct_unclassified", "pct_classified", "n_species"])
        by_comp = {}
        for s in samples:
            u = unclass.get(s)
            classified = round(100.0 - u, 4) if u is not None else ""
            c = comp.get(sid(s), "NA")
            w.writerow([sid(s), c, u if u is not None else "",
                        classified, nsp[s]])
            if u is not None:
                by_comp.setdefault(c, []).append(100.0 - u)

    print("MetaPhlAn %classified (100 - UNCLASSIFIED) by compartment:")
    print(f"{'comp':>16} {'n':>4} {'median':>8}")
    for c in sorted(by_comp):
        v = by_comp[c]
        print(f"{c:>16} {len(v):>4} {st.median(v):>8.2f}")

if __name__ == "__main__":
    main()
