#!/usr/bin/env python3
"""KO x sample abundance table.

Joins the eggNOG KO annotation (catalog_ko_ec.tsv.gz) to the CoverM abundance table
(gene_abundance_long.tsv.gz) and sums raw read counts per KEGG Orthology (bare K#####
ids). A gene annotated with several KOs contributes its full count to each.

Inputs (under $SGMA2_BASE/assembly/genecat): eggnog/catalog_ko_ec.tsv.gz,
gene_abundance_long.tsv.gz, coverm_samples.txt.
Outputs (in eggnog/): ko_abundance_long.tsv.gz, ko_abundance_wide.tsv.gz,
ko_table_report.txt.

MIN_COVERED_FRACTION (environment variable, default 0) drops gene x sample cells
below that covered fraction before summing; outputs then take a _cfNN suffix.
"""
import gzip
import os, sys
from collections import defaultdict

BASE = os.environ["SGMA2_BASE"]
GEN  = f"{BASE}/assembly/genecat"
EGG  = f"{GEN}/eggnog"
KOEC = f"{EGG}/catalog_ko_ec.tsv.gz"
ABUN = f"{GEN}/gene_abundance_long.tsv.gz"
SAMP = f"{GEN}/coverm_samples.txt"
LONG = f"{EGG}/ko_abundance_long.tsv.gz"
WIDE = f"{EGG}/ko_abundance_wide.tsv.gz"
REP  = f"{EGG}/ko_table_report.txt"

MIN_CF = float(os.environ.get("MIN_COVERED_FRACTION", "0"))
if MIN_CF > 0:
    _sfx = "_cf%02d" % int(round(MIN_CF * 100))
    LONG = LONG.replace(".tsv.gz", _sfx + ".tsv.gz")
    WIDE = WIDE.replace(".tsv.gz", _sfx + ".tsv.gz")
    REP  = REP.replace(".txt", _sfx + ".txt")


def norm_ko(field):
    """Split a KEGG_ko field into bare K##### ids; '-' / '' -> none."""
    if not field or field == "-":
        return ()
    out = []
    for tok in field.split(","):
        tok = tok.strip()
        if tok.startswith("ko:"):
            tok = tok[3:]
        if tok and tok != "-":
            out.append(tok)
    return tuple(out)


def main():
    # ---- sample order ----
    with open(SAMP) as fh:
        samples = [ln.strip() for ln in fh if ln.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    print(f"samples: {len(samples)}", flush=True)

    # ---- gene_id -> tuple(KO) for genes that carry >=1 KO ----
    gene_ko = {}
    n_koec = 0
    with gzip.open(KOEC, "rt") as fh:
        header = next(fh).rstrip("\n").split("\t")
        # tolerate either an explicit header or data on line 1
        if header[0] != "gene_id":
            fh.seek(0)
        for ln in fh:
            n_koec += 1
            p = ln.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            kos = norm_ko(p[1])
            if kos:
                gene_ko[p[0]] = kos
    print(f"KO/EC rows: {n_koec} | genes with >=1 KO: {len(gene_ko)}", flush=True)

    # ---- stream abundance, accumulate (ko, sample_idx) -> count ----
    acc = defaultdict(float)
    rows = 0
    used_rows = 0            # abundance rows whose gene carries a KO
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
            kos = gene_ko.get(p[0])
            if kos is None:
                continue
            si = sidx.get(p[1])
            if si is None:
                missing_samples.add(p[1])
                continue
            cnt = float(p[2])
            if cnt == 0:
                continue
            used_rows += 1
            for ko in kos:
                acc[(ko, si)] += cnt
            if rows % 50_000_000 == 0:
                print(f"  ...{rows:,} abundance rows scanned", flush=True)
    print(f"abundance rows: {rows:,} | rows with KO: {used_rows:,}", flush=True)
    if missing_samples:
        print(f"WARNING: {len(missing_samples)} sample ids in abundance not in "
              f"coverm_samples.txt: {sorted(missing_samples)[:10]}", flush=True)

    kos_sorted = sorted({ko for (ko, _) in acc})
    nnz = len(acc)
    print(f"distinct KOs: {len(kos_sorted)} | non-zero (KO,sample) cells: {nnz:,}",
          flush=True)

    def fmt(v):
        return str(int(v)) if v == int(v) else repr(v)

    # ---- long output ----
    with gzip.open(LONG, "wt") as out:
        out.write("ko\tsample_id\tcount\n")
        for ko in kos_sorted:
            for si, s in enumerate(samples):
                v = acc.get((ko, si))
                if v:
                    out.write(f"{ko}\t{s}\t{fmt(v)}\n")

    # ---- wide output (HUMAnN-style: KO rows x sample cols, full sample set) ----
    with gzip.open(WIDE, "wt") as out:
        out.write("ko\t" + "\t".join(samples) + "\n")
        for ko in kos_sorted:
            vals = (acc.get((ko, si), 0) for si in range(len(samples)))
            out.write(ko + "\t" + "\t".join(fmt(v) for v in vals) + "\n")

    # ---- report ----
    with open(REP, "w") as out:
        out.write("SGMA2 KO x sample table report\n")
        out.write(f"samples:                 {len(samples)}\n")
        out.write(f"genes with >=1 KO:       {len(gene_ko)}\n")
        out.write(f"abundance rows scanned:  {rows}\n")
        out.write(f"abundance rows with KO:  {used_rows}\n")
        out.write(f"distinct KOs:            {len(kos_sorted)}\n")
        out.write(f"non-zero KOxsample cells:{nnz}\n")
        out.write(f"multi-KO policy:         full count to each KO\n")
        out.write(f"count filter:            none (raw counts; threshold at analysis)\n")
    print("done", flush=True)


if __name__ == "__main__":
    main()
