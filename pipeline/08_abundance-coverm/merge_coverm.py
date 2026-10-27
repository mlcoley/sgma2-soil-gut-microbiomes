#!/usr/bin/env python3
"""Stack per-sample CoverM outputs into one sparse long gene x sample table.

Inputs (under $SGMA2_BASE/assembly): coverm/<sample>.tsv (Contig, read count,
covered fraction), genecat/catalog_gene_lengths.tsv (gene order reference),
assembly_manifest.tsv.

Outputs: genecat/gene_abundance_long.tsv.gz written to stdout (gene_id, sample_id,
count, covered_fraction; count > 0 rows only), genecat/coverm_samples.txt and
genecat/coverm_merge_report.txt.

Every sample file is checked line by line against the catalog gene order; any
mismatch aborts. Length normalisation is done downstream by joining gene_id to
catalog_gene_lengths.tsv.
"""
import sys, os

BASE     = os.environ["SGMA2_BASE"]
WORK     = os.path.join(BASE, "assembly")
COVERM   = os.path.join(WORK, "coverm")
GENECAT  = os.path.join(WORK, "genecat")
MANIFEST = os.path.join(WORK, "assembly_manifest.tsv")
LENGTHS  = os.path.join(GENECAT, "catalog_gene_lengths.tsv")
SAMPLES_OUT = os.path.join(GENECAT, "coverm_samples.txt")
REPORT_OUT  = os.path.join(GENECAT, "coverm_merge_report.txt")

def log(*a):
    print(*a, file=sys.stderr, flush=True)

# --- gene order reference (col1 of catalog_gene_lengths.tsv, in catalog order) ---
log("Loading gene order from", LENGTHS, "...")
genes = []
with open(LENGTHS) as fh:
    for ln in fh:
        genes.append(ln.split('\t', 1)[0])
NG = len(genes)
log("gene order loaded:", NG, "genes")

# --- samples from the manifest, in manifest order ---
samples = []
with open(MANIFEST) as fh:
    for ln in fh:
        f = ln.rstrip('\n').split('\t')
        if f and f[0]:
            samples.append(f[0])
log("samples in manifest:", len(samples))

# --- pre-flight: every sample TSV must exist before we write anything ---
missing = [s for s in samples if not os.path.isfile(os.path.join(COVERM, s + ".tsv"))]
if missing:
    log("ERROR: missing CoverM TSVs for", len(missing), "samples:", ",".join(missing))
    log("CoverM array incomplete -- aborting before emitting output.")
    sys.exit(1)

out = sys.stdout
report = []
total_nnz = 0

for si, s in enumerate(samples, 1):
    tsv = os.path.join(COVERM, s + ".tsv")
    nnz = 0
    nrow = 0
    with open(tsv) as fh:
        fh.readline()  # skip header
        for i, ln in enumerate(fh):
            g, rest = ln.split('\t', 1)
            if g != genes[i]:
                log("!! ORDER MISMATCH in", s, "at data line", i + 1,
                    "-- expected", genes[i], "got", g)
                sys.exit(2)
            c, cov = rest.split('\t', 1)   # cov still carries its trailing newline
            if c != '0' and float(c) > 0:
                out.write(g); out.write('\t'); out.write(s); out.write('\t')
                out.write(c); out.write('\t'); out.write(cov)
                nnz += 1
            nrow += 1
    if nrow != NG:
        log("!! ROW COUNT MISMATCH in", s, ":", nrow, "data rows vs expected", NG)
        sys.exit(3)
    total_nnz += nnz
    report.append((s, nnz))
    log(f"[{si}/{len(samples)}] {s}: {nnz} nonzero genes ({100.0*nnz/NG:.1f}% of catalog)")

out.flush()

with open(SAMPLES_OUT, "w") as fh:
    for s in samples:
        fh.write(s + "\n")

with open(REPORT_OUT, "w") as fh:
    fh.write("# CoverM merge report\n")
    fh.write(f"# catalog genes: {NG}\n")
    fh.write(f"# samples merged: {len(samples)}\n")
    fh.write(f"# total non-zero cells: {total_nnz}\n")
    fh.write("sample\tnonzero_genes\tpct_of_catalog\n")
    for s, nnz in report:
        fh.write(f"{s}\t{nnz}\t{100.0*nnz/NG:.2f}\n")

log("DONE. total nonzero cells:", total_nnz, "across", len(samples), "samples")
