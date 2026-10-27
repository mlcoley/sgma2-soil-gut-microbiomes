#!/bin/bash
#SBATCH --job-name=arg_gene_abun
#SBATCH --output=logs/arg_gene_abun_%j.out
#SBATCH --error=logs/arg_gene_abun_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Extracts the ARG-bearing catalog genes (gene ids in arg_hits_long.tsv.gz) from the
# gene x sample abundance table, keeping count > 0 rows with covered_fraction.
# Streams the full table, so it runs as a batch job.
# Output: genecat/amr/arg_gene_abundance.tsv.gz (gene_id, sample_id, count,
#   covered_fraction), arg_gene_abundance_report.txt

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
echo "start $(date +'%F %H:%M:%S') on $(hostname)"

BASE="${SGMA2_BASE}"
GEN="${BASE}/assembly/genecat"
HITS="${GEN}/amr/arg_hits_long.tsv.gz"
ABUN="${GEN}/gene_abundance_long.tsv.gz"
OUT="${GEN}/amr/arg_gene_abundance.tsv.gz"
REP="${GEN}/amr/arg_gene_abundance_report.txt"
mkdir -p logs "${GEN}/amr"

[ -f "$HITS" ] || { echo "ERROR: missing $HITS"; exit 1; }
[ -f "$ABUN" ] || { echo "ERROR: missing $ABUN"; exit 1; }

TMP="${GEN}/amr/.arg_gene_abun.${SLURM_JOB_ID}"
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT

# ---- the ARG-bearing catalog gene ids -------------------------------
zcat "$HITS" | awk -F'\t' 'NR>1 && $1!=""{print $1}' | sort -u > "${TMP}/arg_gene_ids.txt"
n_ids=$(wc -l < "${TMP}/arg_gene_ids.txt")
echo "ARG-bearing catalog genes: ${n_ids}"

# ---- stream the abundance matrix, keep ARG genes with count>0 -------------
# gene_abundance_long may or may not carry a header row (build_arg_table.py
# tolerates both), so detect rather than assume: a first field literally equal
# to "gene_id" is a header.
zcat "$ABUN" | awk -F'\t' -v OFS='\t' -v IDS="${TMP}/arg_gene_ids.txt" '
  BEGIN{ while ((getline g < IDS) > 0) keep[g]=1
         print "gene_id","sample_id","count","covered_fraction" }
  FNR==1 && $1=="gene_id" { next }
  ($1 in keep) && ($3+0) > 0 { print $1,$2,$3,($4==""?"NA":$4); n++ }
  END{ print "rows_written="n > "/dev/stderr" }
' | gzip > "${TMP}/out.tsv.gz"

# ---- report + atomic publish ----------------------------------------------
rows=$(zcat "${TMP}/out.tsv.gz" | tail -n +2 | wc -l)
nsamp=$(zcat "${TMP}/out.tsv.gz" | tail -n +2 | cut -f2 | sort -u | wc -l)
ngene=$(zcat "${TMP}/out.tsv.gz" | tail -n +2 | cut -f1 | sort -u | wc -l)
{
  echo "ARG-gene abundance slice — $(date +'%F %H:%M:%S')"
  echo "ARG-bearing catalog genes (from arg_hits_long): ${n_ids}"
  echo "non-zero gene x sample rows written:            ${rows}"
  echo "distinct samples represented:                   ${nsamp}"
  echo "distinct ARG genes with >=1 read somewhere:     ${ngene}"
} | tee "$REP"

mv "${TMP}/out.tsv.gz" "${OUT}.partial"
mv "${OUT}.partial" "$OUT"
echo "wrote $OUT"
echo "done $(date +'%F %H:%M:%S')"
