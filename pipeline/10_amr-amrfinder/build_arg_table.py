#!/usr/bin/env python3
"""ARG x sample abundance table.

Joins the harmonised AMR hit table (arg_hits_long.tsv.gz) to the CoverM abundance
table (gene_abundance_long.tsv.gz) and sums raw read counts per ARG gene symbol. A gene
with several ARG hits contributes its full count to each (tool, arg_symbol). No scope
filter is applied; arg_annotation.tsv carries element type, scope, drug class and
subclass per symbol for filtering downstream.

Inputs (under $SGMA2_BASE/assembly/genecat): amr/arg_hits_long.tsv.gz,
gene_abundance_long.tsv.gz, coverm_samples.txt.
Outputs (in amr/): arg_abundance_long.tsv.gz, arg_abundance_wide.tsv.gz,
arg_annotation.tsv, arg_table_report.txt.

MIN_COVERED_FRACTION (environment variable, default 0) drops gene x sample cells
below that covered fraction before summing; outputs then take a _cfNN suffix.
"""
import gzip
import os, sys
from collections import defaultdict

BASE = os.environ["SGMA2_BASE"]
GEN  = f"{BASE}/assembly/genecat"
AMR  = f"{GEN}/amr"
HITS = f"{AMR}/arg_hits_long.tsv.gz"
ABUN = f"{GEN}/gene_abundance_long.tsv.gz"
SAMP = f"{GEN}/coverm_samples.txt"
LONG = f"{AMR}/arg_abundance_long.tsv.gz"
WIDE = f"{AMR}/arg_abundance_wide.tsv.gz"
ANNO = f"{AMR}/arg_annotation.tsv"
REP  = f"{AMR}/arg_table_report.txt"

MIN_CF = float(os.environ.get("MIN_COVERED_FRACTION", "0"))
if MIN_CF > 0:
    _sfx = "_cf%02d" % int(round(MIN_CF * 100))
    LONG = LONG.replace(".tsv.gz", _sfx + ".tsv.gz")
    WIDE = WIDE.replace(".tsv.gz", _sfx + ".tsv.gz")
    REP  = REP.replace(".txt", _sfx + ".txt")



def main():
    # ---- sample order ----
    with open(SAMP) as fh:
        samples = [ln.strip() for ln in fh if ln.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    print(f"samples: {len(samples)}", flush=True)

    # ---- gene_id -> set((tool, arg_symbol)); per-symbol annotation ----
    # arg_hits_long header:
    #   gene_id tool arg_symbol arg_name scope element_type element_subtype
    #   drug_class subclass method pct_coverage pct_identity ref_accession
    gene_args = defaultdict(set)
    anno = {}                      # arg_symbol -> dict of metadata (first-seen)
    n_hits = 0
    with gzip.open(HITS, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        need = ["gene_id", "tool", "arg_symbol"]
        for c in need:
            if c not in col:
                sys.exit(f"ERROR: arg_hits_long missing column '{c}'; header={header}")
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) <= col["arg_symbol"]:
                continue
            gid = p[col["gene_id"]]
            tool = p[col["tool"]]
            sym = p[col["arg_symbol"]]
            if not sym:
                continue
            n_hits += 1
            gene_args[gid].add((tool, sym))
            if sym not in anno:
                anno[sym] = {
                    "tool":            tool,
                    "element_type":    p[col["element_type"]]    if "element_type"    in col and col["element_type"]    < len(p) else "",
                    "scope":           p[col["scope"]]           if "scope"           in col and col["scope"]           < len(p) else "",
                    "drug_class":      p[col["drug_class"]]      if "drug_class"      in col and col["drug_class"]      < len(p) else "",
                    "subclass":        p[col["subclass"]]        if "subclass"        in col and col["subclass"]        < len(p) else "",
                    "n_genes":         0,
                }
    # count distinct genes per symbol (for the sidecar / QC)
    for gid, pairs in gene_args.items():
        for _tool, sym in pairs:
            anno[sym]["n_genes"] += 1
    print(f"ARG hits: {n_hits} | genes with >=1 ARG: {len(gene_args)} | "
          f"distinct symbols: {len(anno)}", flush=True)

    # ---- stream abundance, accumulate (tool, arg_symbol, sample_idx) -> count ----
    acc = defaultdict(float)
    rows = 0
    used_rows = 0
    missing_samples = set()
    with gzip.open(ABUN, "rt") as fh:
        first = fh.readline().rstrip("\n").split("\t")
        if MIN_CF > 0 and first and first[0] == "gene_id" and len(first) < 4:
            raise SystemExit("MIN_COVERED_FRACTION set but %s lacks a covered_fraction column" % ABUN)
        if first and first[0] != "gene_id":       # line 1 was data, not a header
            fh.seek(0)
        for ln in fh:
            rows += 1
            p = ln.rstrip("\n").split("\t")
            if MIN_CF > 0 and float(p[3]) < MIN_CF:
                continue          # breadth filter: skip low-coverage gene x sample cells
            if len(p) < 3:
                continue
            pairs = gene_args.get(p[0])
            if pairs is None:
                continue
            si = sidx.get(p[1])
            if si is None:
                missing_samples.add(p[1])
                continue
            cnt = float(p[2])
            if cnt == 0:
                continue
            used_rows += 1
            for tool, sym in pairs:
                acc[(tool, sym, si)] += cnt
            if rows % 50_000_000 == 0:
                print(f"  ...{rows:,} abundance rows scanned", flush=True)
    print(f"abundance rows: {rows:,} | rows with ARG: {used_rows:,}", flush=True)
    if missing_samples:
        print(f"WARNING: {len(missing_samples)} sample ids in abundance not in "
              f"coverm_samples.txt: {sorted(missing_samples)[:10]}", flush=True)

    syms_sorted = sorted({sym for (_t, sym, _s) in acc})
    tools = sorted({t for (t, _s, _i) in acc})
    nnz = len(acc)
    print(f"distinct ARG symbols (non-zero): {len(syms_sorted)} | tools: {tools} | "
          f"non-zero cells: {nnz:,}", flush=True)

    def fmt(v):
        return str(int(v)) if v == int(v) else repr(v)

    # ---- long output (tool-tagged) ----
    with gzip.open(LONG, "wt") as out:
        out.write("tool\targ_symbol\tsample_id\tcount\n")
        for (tool, sym, si), v in sorted(acc.items()):
            if v:
                out.write(f"{tool}\t{sym}\t{samples[si]}\t{fmt(v)}\n")

    # ---- wide output (arg_symbol x sample). Single tool now; if a future tool
    #      re-uses a symbol, counts for that symbol are summed across tools here
    #      (kept separable in the long table). Symbol nomenclatures differ across
    #      tools so collisions are near-zero; flagged in the report if any. ----
    sym_tools = defaultdict(set)
    for (tool, sym, _si) in acc:
        sym_tools[sym].add(tool)
    collisions = [s for s, t in sym_tools.items() if len(t) > 1]
    wide_acc = defaultdict(float)
    for (tool, sym, si), v in acc.items():
        wide_acc[(sym, si)] += v
    with gzip.open(WIDE, "wt") as out:
        out.write("arg_symbol\t" + "\t".join(samples) + "\n")
        for sym in syms_sorted:
            vals = (wide_acc.get((sym, si), 0) for si in range(len(samples)))
            out.write(sym + "\t" + "\t".join(fmt(v) for v in vals) + "\n")

    # ---- per-symbol annotation sidecar (drives analysis-layer filtering) ----
    with open(ANNO, "w") as out:
        out.write("arg_symbol\ttool\telement_type\tscope\tdrug_class\tsubclass\tn_genes\n")
        for sym in syms_sorted:
            a = anno.get(sym, {})
            out.write("\t".join([
                sym, a.get("tool", ""), a.get("element_type", ""),
                a.get("scope", ""), a.get("drug_class", ""),
                a.get("subclass", ""), str(a.get("n_genes", 0)),
            ]) + "\n")

    # ---- report ----
    with open(REP, "w") as out:
        out.write("SGMA2 ARG x sample table report\n")
        out.write(f"samples:                  {len(samples)}\n")
        out.write(f"tools:                    {tools}\n")
        out.write(f"ARG hits (arg_hits_long): {n_hits}\n")
        out.write(f"genes with >=1 ARG:       {len(gene_args)}\n")
        out.write(f"abundance rows scanned:   {rows}\n")
        out.write(f"abundance rows with ARG:  {used_rows}\n")
        out.write(f"distinct ARG symbols:     {len(anno)}\n")
        out.write(f"  ...with non-zero abund: {len(syms_sorted)}\n")
        out.write(f"non-zero ARGxsample cells:{nnz}\n")
        out.write(f"multi-hit policy:         full count to each ARG symbol\n")
        out.write(f"count filter:             none (raw counts; threshold at analysis)\n")
        out.write(f"scope filter:             none (element_type/scope in arg_annotation.tsv)\n")
        if collisions:
            out.write(f"\nWARNING: {len(collisions)} symbols called by >1 tool "
                      f"(summed in wide, separable in long): {collisions[:20]}\n")
    print("done", flush=True)


if __name__ == "__main__":
    main()
