#!/usr/bin/env python3
"""Merge per-chunk AMRFinderPlus tables.

Writes a native concatenated table and a harmonised long table keyed on gene_id with
a `tool` column, so another AMR tool can be appended with the same schema. Columns
are mapped by header name, not position.

Inputs: genecat/amr/amrfinder/chunk_*.tsv (under $SGMA2_BASE/assembly).
Outputs: genecat/amr/amrfinder/catalog_amrfinder.tsv.gz,
genecat/amr/arg_hits_long.tsv.gz, genecat/amr/amr_merge_report.txt.
"""
import gzip, glob, os, sys
from collections import Counter

BASE = os.environ["SGMA2_BASE"]
AMR  = f"{BASE}/assembly/genecat/amr"
INDIR = f"{AMR}/amrfinder"
FULL  = f"{AMR}/amrfinder/catalog_amrfinder.tsv.gz"
LONG  = f"{AMR}/arg_hits_long.tsv.gz"
REP   = f"{AMR}/amr_merge_report.txt"

# harmonized_col -> candidate AMRFinderPlus header names (version-tolerant).
# AMRFinderPlus v4 renamed several columns from v3; both spellings are listed.
MAP = {
    "gene_id":         ["Protein identifier", "Protein id"],
    "arg_symbol":      ["Element symbol", "Gene symbol"],
    "arg_name":        ["Element name", "Sequence name"],
    "scope":           ["Scope"],
    "element_type":    ["Type", "Element type"],
    "element_subtype": ["Subtype", "Element subtype"],
    "drug_class":      ["Class"],
    "subclass":        ["Subclass"],
    "method":          ["Method"],
    "pct_coverage":    ["% Coverage of reference", "% Coverage of reference sequence"],
    "pct_identity":    ["% Identity to reference", "% Identity to reference sequence"],
    "ref_accession":   ["Closest reference accession", "Accession of closest sequence"],
}
HARM_COLS = ["gene_id", "tool"] + [c for c in MAP if c != "gene_id"]


def resolve(header):
    """Map harmonized column -> index in this file's header, by name."""
    idx = {}
    for hc, cands in MAP.items():
        idx[hc] = next((header.index(c) for c in cands if c in header), None)
    return idx


def main():
    files = sorted(glob.glob(f"{INDIR}/chunk_*.tsv"))
    if not files:
        sys.exit(f"ERROR: no chunk_*.tsv in {INDIR}")
    print(f"chunk files: {len(files)}", flush=True)

    full = gzip.open(FULL, "wt")
    lng  = gzip.open(LONG, "wt")
    lng.write("\t".join(HARM_COLS) + "\n")

    full_header = None
    total = 0
    per_class = Counter()
    per_symbol = Counter()
    genes_with_arg = set()

    for i, fp in enumerate(files):
        with open(fp) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            if full_header is None:                 # keep one native header
                full_header = header
                full.write("\t".join(header) + "\n")
            idx = resolve(header)
            if idx["gene_id"] is None:
                sys.exit(f"ERROR: no gene_id column in {fp} header: {header}")
            for ln in fh:
                row = ln.rstrip("\n").split("\t")
                full.write(ln if ln.endswith("\n") else ln + "\n")
                total += 1
                gid = row[idx["gene_id"]]
                genes_with_arg.add(gid)
                out = [gid, "amrfinder"]
                for hc in HARM_COLS[2:]:
                    j = idx[hc]
                    out.append(row[j] if (j is not None and j < len(row)) else "")
                lng.write("\t".join(out) + "\n")
                if idx["drug_class"] is not None and idx["drug_class"] < len(row):
                    per_class[row[idx["drug_class"]]] += 1
                if idx["arg_symbol"] is not None and idx["arg_symbol"] < len(row):
                    per_symbol[row[idx["arg_symbol"]]] += 1

    full.close(); lng.close()

    with open(REP, "w") as out:
        out.write("SGMA2 AMR merge report (AMRFinderPlus)\n")
        out.write(f"chunks merged:        {len(files)}\n")
        out.write(f"total ARG hits:       {total}\n")
        out.write(f"distinct genes w/ ARG:{len(genes_with_arg)}\n")
        out.write(f"distinct gene symbols:{len(per_symbol)}\n")
        out.write("\ntop 20 drug classes (by hit count):\n")
        for cls, n in per_class.most_common(20):
            out.write(f"  {n:>8}  {cls}\n")
        out.write("\ntop 20 gene symbols (by hit count):\n")
        for sym, n in per_symbol.most_common(20):
            out.write(f"  {n:>8}  {sym}\n")
    print(f"done: {total} hits, {len(genes_with_arg)} genes with ARG", flush=True)


if __name__ == "__main__":
    main()
