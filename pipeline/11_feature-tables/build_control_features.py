#!/usr/bin/env python3
"""Control libraries in ARG-symbol and KO feature space.

Turns the per-control CoverM outputs into feature profiles. Writes an audit table
(every feature x control with read counts), a proposed exclusion list from negative and
PCR blanks only (positive duplicates and mocks never contribute), and an impact table
showing how many features each read threshold would remove. The proposed exclusion list
lists features at >= 5 reads (DEFAULT_MIN_READS) and is not applied until copied to
feature_tables/ without ".unapplied" in its name.

Inputs (under $SGMA2_BASE/assembly): coverm/<control>.tsv, control_manifest.tsv
(sample, control_type), genecat/amr/arg_hits_long.tsv.gz,
genecat/eggnog/catalog_ko_ec.tsv.gz.
Outputs (genecat/controls/): control_features_detected.tsv,
control_features_exclusion_list.unapplied.tsv, control_exclusion_list_impact.tsv,
control_features_report.txt.
"""
import gzip, os, sys
from collections import defaultdict

BASE = os.environ["SGMA2_BASE"]
WORK = f"{BASE}/assembly"
GEN  = f"{WORK}/genecat"
COV  = f"{WORK}/coverm"
MAN  = f"{WORK}/control_manifest.tsv"
HITS = f"{GEN}/amr/arg_hits_long.tsv.gz"
KOEC = f"{GEN}/eggnog/catalog_ko_ec.tsv.gz"
OUTD = f"{GEN}/controls"

# Detection threshold (reads) for the proposed list; the impact table spans a range.
DEFAULT_MIN_READS = 5
IMPACT_THRESHOLDS = [1, 2, 5, 10, 50]
# Only these control types contribute to an exclusion list.
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


def load_gene_to_ko():
    # catalog_ko_ec.tsv.gz: gene_id, KEGG_ko, EC. KEGG_ko is a comma list with a
    # 'ko:' prefix; '-' means none. Same parsing as build_ko_table.py
    g2k = {}
    with gzip.open(KOEC, "rt") as fh:
        first = fh.readline().rstrip("\n").split("\t")
        if first and first[0] not in ("gene_id", "#query"):
            fh.seek(0)
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            raw = p[1]
            if not raw or raw == "-":
                continue
            kos = tuple(k.split(":")[-1] for k in raw.split(",") if k and k != "-")
            if kos:
                g2k[p[0]] = kos
    return g2k


def main():
    os.makedirs(OUTD, exist_ok=True)
    controls = load_manifest()
    print(f"controls in manifest: {len(controls)}", flush=True)

    present = [(s, t) for s, t in controls if os.path.exists(f"{COV}/{s}.tsv")]
    missing = [s for s, _ in controls if not os.path.exists(f"{COV}/{s}.tsv")]
    if missing:
        print(f"WARNING: no CoverM output for {len(missing)}: {','.join(missing)}", flush=True)
    if not present:
        sys.exit("ERROR: no control CoverM outputs found - run coverm_array.sh with control_manifest.tsv first")

    print("loading gene -> ARG ...", flush=True)
    g2a = load_gene_to_arg()
    print(f"  ARG-bearing genes: {len(g2a)}", flush=True)
    print("loading gene -> KO (large) ...", flush=True)
    g2k = load_gene_to_ko()
    print(f"  KO-annotated genes: {len(g2k)}", flush=True)

    # feature -> control -> summed count
    detected = defaultdict(lambda: defaultdict(int))
    ctype = dict(present)
    for s, t in present:
        n_nz = 0
        with open(f"{COV}/{s}.tsv") as fh:
            fh.readline()                       # header
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                if len(p) < 2:
                    continue
                try:
                    c = float(p[1])
                except ValueError:
                    continue
                if c <= 0:
                    continue
                n_nz += 1
                gid = p[0]
                for sym in g2a.get(gid, ()):
                    detected[("ARG", sym)][s] += int(c)
                for ko in g2k.get(gid, ()):
                    detected[("KO", ko)][s] += int(c)
        print(f"  {s} ({t}): {n_nz} genes with >=1 read", flush=True)

    # ---- audit table -------------------------------------------------------
    aud = f"{OUTD}/control_features_detected.tsv"
    with open(aud, "w") as out:
        out.write("feature_type\tfeature\tcontrol\tcontrol_type\tcount\n")
        for (ft, feat), per in sorted(detected.items()):
            for s, c in sorted(per.items()):
                out.write(f"{ft}\t{feat}\t{s}\t{ctype[s]}\t{c}\n")
    print(f"wrote {aud}", flush=True)

    # ---- exclusion list impact across thresholds --------------------------------
    blk_controls = {s for s, t in present if t in EXCLUSION_LIST_TYPES}
    imp = f"{OUTD}/control_exclusion_list_impact.tsv"
    with open(imp, "w") as out:
        out.write("min_reads\tfeature_type\tn_features_excluded\n")
        for thr in IMPACT_THRESHOLDS:
            cnt = defaultdict(int)
            for (ft, feat), per in detected.items():
                if any(c >= thr for s, c in per.items() if s in blk_controls):
                    cnt[ft] += 1
            for ft in ("ARG", "KO"):
                out.write(f"{thr}\t{ft}\t{cnt[ft]}\n")
    print(f"wrote {imp}", flush=True)

    # ---- unapplied exclusion list (requires a deliberate rename to activate) ----
    cand = f"{OUTD}/control_features_exclusion_list.unapplied.tsv"
    n_c = 0
    with open(cand, "w") as out:
        out.write("feature\tfeature_type\tmax_count_in_negative_control\tn_negative_controls\n")
        for (ft, feat), per in sorted(detected.items()):
            hits = {s: c for s, c in per.items()
                    if s in blk_controls and c >= DEFAULT_MIN_READS}
            if hits:
                out.write(f"{feat}\t{ft}\t{max(hits.values())}\t{len(hits)}\n")
                n_c += 1
    print(f"wrote {cand} ({n_c} features proposed for exclusion at >={DEFAULT_MIN_READS} read)", flush=True)

    with open(f"{OUTD}/control_features_report.txt", "w") as out:
        out.write("Control feature profiles\n\n")
        out.write(f"controls mapped: {len(present)} of {len(controls)}\n")
        if missing:
            out.write(f"missing CoverM output: {','.join(missing)}\n")
        out.write(f"negative controls used for the proposed exclusion list: "
                  f"{','.join(sorted(blk_controls))}\n")
        out.write(f"features proposed for exclusion at >={DEFAULT_MIN_READS} read: {n_c}\n\n")
        out.write("Positive duplicates and mock are NOT used for contamination filtering - they are real\n"
                  "material, carried here for FDR calibration and recovery checks.\n\n")
        out.write("To apply the contamination filter, choose a threshold from control_exclusion_list_impact.tsv and\n")
        out.write("write the chosen features, one per line under a `feature` header, to\n")
        out.write("feature_tables/control_features_exclusion_list.tsv\n")
    print("done", flush=True)


if __name__ == "__main__":
    main()
