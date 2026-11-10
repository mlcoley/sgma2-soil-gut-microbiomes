#!/usr/bin/env python3
"""Control detections at catalog-gene level, restricted to ARG-bearing genes.

Same design as build_control_features.py, keyed on gene_id instead of ARG symbol:
only negative and PCR blanks contribute to the proposed exclusion list (genes at >= 5
reads), which is not applied until copied to feature_tables/ without ".unapplied" in its name.

Inputs (under $SGMA2_BASE/assembly): coverm/<control>.tsv, control_manifest.tsv
(sample, control_type), genecat/amr/arg_hits_long.tsv.gz.
Outputs (genecat/controls/): control_arg_genes_detected.tsv,
control_gene_exclusion_list_impact.tsv, control_arg_genes_exclusion_list.unapplied.tsv,
control_gene_features_report.txt.
"""
import gzip, os, sys
from collections import defaultdict

BASE = os.environ["SGMA2_BASE"]
WORK = f"{BASE}/assembly"
GEN  = f"{WORK}/genecat"
COV  = f"{WORK}/coverm"
MAN  = f"{WORK}/control_manifest.tsv"
HITS = f"{GEN}/amr/arg_hits_long.tsv.gz"
OUTD = f"{GEN}/controls"

IMPACT_THRESHOLDS = [1, 2, 5, 10, 50]
DEFAULT_MIN_READS = 5
EXCLUSION_LIST_TYPES = {"neg_blank", "pcr_blank"}


def load_manifest():
    rows = []
    with open(MAN) as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0]:
                rows.append((p[0], p[1]))
    return rows


def load_gene_to_arg():
    g2a = defaultdict(set)
    with gzip.open(HITS, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) > col["arg_symbol"] and p[col["arg_symbol"]]:
                g2a[p[col["gene_id"]]].add(p[col["arg_symbol"]])
    return g2a


def main():
    os.makedirs(OUTD, exist_ok=True)
    controls = load_manifest()
    present = [(s, t) for s, t in controls if os.path.exists(f"{COV}/{s}.tsv")]
    missing = [s for s, _ in controls if not os.path.exists(f"{COV}/{s}.tsv")]
    if missing:
        print(f"WARNING: no CoverM output for {len(missing)}: {','.join(missing)}", flush=True)
    if not present:
        sys.exit("ERROR: no control CoverM outputs found - run coverm_array.sh with control_manifest.tsv first")
    print(f"controls with CoverM output: {len(present)} of {len(controls)}", flush=True)

    g2a = load_gene_to_arg()
    print(f"ARG-bearing catalog genes: {len(g2a)}", flush=True)

    ctype = dict(present)
    detected = defaultdict(lambda: defaultdict(int))   # gene_id -> control -> count
    for s, t in present:
        n_hit = 0
        with open(f"{COV}/{s}.tsv") as fh:
            fh.readline()                              # header
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                if len(p) < 2:
                    continue
                gid = p[0]
                if gid not in g2a:
                    continue
                try:
                    c = float(p[1])
                except ValueError:
                    continue
                if c <= 0:
                    continue
                detected[gid][s] += int(c)
                n_hit += 1
        print(f"  {s} ({t}): {n_hit} ARG genes with >=1 read", flush=True)

    # ---- audit table -------------------------------------------------------
    aud = f"{OUTD}/control_arg_genes_detected.tsv"
    with open(aud, "w") as out:
        out.write("gene_id\targ_symbol\tcontrol\tcontrol_type\tcount\n")
        for gid, per in sorted(detected.items()):
            syms = ",".join(sorted(g2a[gid]))
            for s, c in sorted(per.items()):
                out.write(f"{gid}\t{syms}\t{s}\t{ctype[s]}\t{c}\n")
    print(f"wrote {aud}", flush=True)

    # ---- impact across thresholds ------------------------------------------
    blk = {s for s, t in present if t in EXCLUSION_LIST_TYPES}
    def excluded_at(thr):
        return {gid for gid, per in detected.items()
                if max((c for s, c in per.items() if s in blk), default=0) >= thr}
    imp = f"{OUTD}/control_gene_exclusion_list_impact.tsv"
    with open(imp, "w") as out:
        out.write("min_reads\tn_genes_excluded\tpct_of_arg_genes\tn_symbols_affected\n")
        for thr in IMPACT_THRESHOLDS:
            g = excluded_at(thr)
            syms = {x for gid in g for x in g2a[gid]}
            out.write(f"{thr}\t{len(g)}\t{100.0 * len(g) / max(len(g2a), 1):.2f}\t{len(syms)}\n")
    print(f"wrote {imp}", flush=True)

    # ---- unapplied exclusion list (deliberate rename required to activate) ------
    cand = f"{OUTD}/control_arg_genes_exclusion_list.unapplied.tsv"
    with open(cand, "w") as out:
        out.write("gene_id\targ_symbol\tmax_count_in_negative_control\tn_negative_controls\n")
        for gid in sorted(excluded_at(DEFAULT_MIN_READS)):
            per = {s: c for s, c in detected[gid].items() if s in blk}
            out.write(f"{gid}\t{','.join(sorted(g2a[gid]))}\t{max(per.values())}\t{len(per)}\n")
    print(f"wrote {cand}", flush=True)

    rep = f"{OUTD}/control_gene_features_report.txt"
    with open(rep, "w") as out:
        out.write("Control detections at catalog-gene level\n\n")
        out.write(f"controls mapped: {len(present)} of {len(controls)}\n")
        out.write(f"ARG-bearing catalog genes: {len(g2a)}\n")
        out.write(f"negative controls used for the proposed exclusion list: "
                  f"{','.join(sorted(blk))}\n")
        out.write(f"genes proposed for exclusion at >={DEFAULT_MIN_READS} reads: {len(excluded_at(DEFAULT_MIN_READS))}\n")
        for thr in IMPACT_THRESHOLDS:
            out.write(f"  genes excluded at >={thr} read(s): {len(excluded_at(thr))}\n")
        out.write("\nPositive duplicates and mock are NOT used for contamination filtering.\n")
        out.write("\nTo apply the gene-level contamination filter, choose a threshold from\n")
        out.write("control_gene_exclusion_list_impact.tsv and write the chosen gene_ids, one per\n")
        out.write("line under a `gene_id` header, to feature_tables/control_arg_genes_exclusion_list.tsv\n")
    print(f"wrote {rep}", flush=True)
    print("done.", flush=True)


if __name__ == "__main__":
    main()
