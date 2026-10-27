#!/usr/bin/env python3
"""Per-gene summary over the gene catalog.

Writes one row per detected catalog gene, restricted to the analysed samples listed
in --cohort: prev_{gut,dung,soil} (samples with count >= 1),
prevcf50_{gut,dung,soil} (samples with covered_fraction >= 0.5), total_count,
mean_nonzero, gene_length, og_narrow (narrowest eggNOG orthologous group), og_bact
(Bacteria-level OG, where present) and ko (first listed KEGG KO).

Genes called by AMRFinderPlus (--exclude) and genes without an eggNOG annotation are
dropped. Matching and semantic exclusions are done downstream.

gene_abundance_long.tsv.gz has no header; a header row is detected, not assumed.
"""

import argparse, gzip, os, sys
from itertools import chain
from sys import intern


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def log(msg):
    print(msg, flush=True)


def read_cohort(path):
    """sample_id -> compartment index (0 gut, 1 dung, 2 soil)."""
    idx = {"gut": 0, "dung": 1, "soil": 2}
    comp = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        assert header[:2] == ["sample_id", "compartment"], f"unexpected header: {header}"
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            comp[p[0]] = idx[p[1]]
    return comp


def read_excludes(path):
    with open(path) as fh:
        return {ln.strip() for ln in fh if ln.strip()}


def parse_ogs(field):
    """eggNOG_OGs is 'OG@taxid|name,OG@taxid|name,...' ordered root -> narrowest.

    Returns (narrowest, bacteria_level). Either may be '-'.
    """
    if not field or field == "-":
        return "-", "-"
    parts = field.split(",")
    narrow = parts[-1].split("@")[0]
    bact = "-"
    for p in parts:
        # match on the taxid, not the display name — names vary by eggNOG build
        if "@2|" in p:
            bact = p.split("@")[0]
            break
    return narrow, bact


def read_eggnog(path, excludes):
    """gene_id -> (og_narrow, og_bact, ko), skipping excluded genes.

    The file is per-chunk emapper .annotations concatenated, so '#' header and
    comment blocks recur throughout. Column positions are read from the first
    '#query' header seen and asserted against the emapper v2 layout.
    """
    ann = {}
    col_ogs = col_ko = None
    n_lines = 0
    with opener(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                if col_ogs is None and ln.lstrip("#").startswith("query"):
                    hdr = ln.lstrip("#").rstrip("\n").split("\t")
                    col_ogs = hdr.index("eggNOG_OGs")
                    col_ko = hdr.index("KEGG_ko")
                    log(f"  eggNOG header found: eggNOG_OGs at {col_ogs}, KEGG_ko at {col_ko}")
                continue
            if col_ogs is None:
                # emapper v2 default layout, used only if no header was present
                col_ogs, col_ko = 4, 11
                log("  WARNING: no '#query' header seen; falling back to emapper v2 columns 4/11")
            p = ln.rstrip("\n").split("\t")
            if len(p) <= col_ko:
                continue
            g = p[0]
            if g in excludes:
                continue
            narrow, bact = parse_ogs(p[col_ogs])
            ko = p[col_ko]
            ko = ko.split(",")[0].replace("ko:", "") if ko and ko != "-" else "-"
            # OG and KO strings repeat across millions of genes — interning
            # collapses them to one object each and saves several GB.
            ann[g] = (intern(narrow), intern(bact), intern(ko))
            n_lines += 1
            if n_lines % 5_000_000 == 0:
                log(f"  eggNOG rows read: {n_lines:,}  retained: {len(ann):,}")
    return ann


def read_lengths(path, keep):
    lens = {}
    with opener(path) as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            if p[0] in keep:
                try:
                    lens[p[0]] = int(p[1])
                except ValueError:
                    continue
    return lens


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--abundance", required=True)
    ap.add_argument("--lengths", required=True)
    ap.add_argument("--eggnog", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--exclude", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-cohort", type=int, default=194)
    a = ap.parse_args()

    log("== per-gene catalog summary ==")

    comp = read_cohort(a.cohort)
    log(f"cohort samples: {len(comp)}")
    assert len(comp) == a.n_cohort, f"expected {a.n_cohort} cohort samples, got {len(comp)}"

    excludes = read_excludes(a.exclude)
    log(f"AMRFinderPlus genes excluded: {len(excludes):,}")

    log("reading eggNOG annotations ...")
    ann = read_eggnog(a.eggnog, excludes)
    log(f"annotated, non-ARG catalog genes: {len(ann):,}")
    assert ann, "no eggNOG annotations parsed — check the header/column detection"

    log("reading gene lengths ...")
    lens = read_lengths(a.lengths, ann)
    log(f"lengths found for {len(lens):,} of {len(ann):,} genes")

    log("streaming gene abundance ...")
    # per gene: [p_gut, p_dung, p_soil, cf_gut, cf_dung, cf_soil, total, n_nonzero]
    acc = {}
    n_rows = n_kept = 0
    with opener(a.abundance) as fh:
        # gene_abundance_long.tsv.gz has no header; detect rather than assume.
        first = fh.readline()
        p0 = first.rstrip("\n").split("\t")
        if p0[:1] == ["gene_id"]:
            assert p0[:4] == ["gene_id", "sample_id", "count", "covered_fraction"], \
                f"header row present but unexpected: {p0}"
            log("  header row detected and consumed")
            stream = fh
        else:
            assert len(p0) >= 4 and p0[2].isdigit(), \
                f"first line is neither a header nor a 4-column data row: {p0}"
            log("  no header row — treating line 1 as data (expected for this file)")
            stream = chain([first], fh)

        for ln in stream:
            n_rows += 1
            if n_rows % 20_000_000 == 0:
                log(f"  abundance rows: {n_rows:,}  genes accumulated: {len(acc):,}")
            p = ln.rstrip("\n").split("\t")
            g = p[0]
            rec = acc.get(g)
            if rec is None:
                if g not in ann:
                    continue
                c = comp.get(p[1])
                if c is None:
                    continue
                rec = acc[g] = [0, 0, 0, 0, 0, 0, 0, 0]
            else:
                c = comp.get(p[1])
                if c is None:
                    continue
            n_kept += 1
            rec[c] += 1
            if float(p[3]) >= 0.5:
                rec[3 + c] += 1
            rec[6] += int(p[2])
            rec[7] += 1
    log(f"abundance rows read {n_rows:,}, retained {n_kept:,}, genes detected {len(acc):,}")

    log(f"writing {a.out} ...")
    tmp = a.out + ".partial"
    n_out = n_nolen = 0
    with gzip.open(tmp, "wt") as out:
        out.write("gene_id\tprev_gut\tprev_dung\tprev_soil\t"
                  "prevcf50_gut\tprevcf50_dung\tprevcf50_soil\t"
                  "total_count\tmean_nonzero\tgene_length\tog_narrow\tog_bact\tko\n")
        for g, r in acc.items():
            L = lens.get(g)
            if L is None:
                n_nolen += 1
                continue
            narrow, bact, ko = ann[g]
            mean_nz = r[6] / r[7] if r[7] else 0.0
            out.write(f"{g}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t"
                      f"{r[6]}\t{mean_nz:.4f}\t{L}\t{narrow}\t{bact}\t{ko}\n")
            n_out += 1
    os.replace(tmp, a.out)
    log(f"wrote {n_out:,} gene rows ({n_nolen:,} dropped for missing length)")

    # A silently-empty output here would look like a clean run and waste a
    # second submission before anyone noticed.
    assert n_out > 0, "no rows written — refusing to publish an empty summary"
    log("done.")


if __name__ == "__main__":
    main()
