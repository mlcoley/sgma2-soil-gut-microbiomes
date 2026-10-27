#!/bin/bash
#SBATCH --job-name=mpa_merge
#SBATCH --output=logs/mpa_merge_%j.out
#SBATCH --error=logs/mpa_merge_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Merges per-sample MetaPhlAn profiles with merge_metaphlan_tables.py, splits
# species and genus tables, and writes the classified summary.
# Run after the MetaPhlAn array finishes.

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
set +u
source ~/.bashrc
conda activate biobakery
set -u

BASE="${SGMA2_BASE}"
MPADIR="${BASE}/biobakery/metaphlan"
OUT="${BASE}/biobakery/merged_metaphlan"
ACC="${BASE}/posthost_qc/sgma2_step1_read_accounting.tsv"
mkdir -p "$OUT" logs

echo "=== merge $(date +'%F %H:%M:%S') ==="
merge_metaphlan_tables.py "${MPADIR}"/*_metaphlan.txt > "${OUT}/metaphlan_merged_all.tsv"
n_samples=$(($(head -1 "${OUT}/metaphlan_merged_all.tsv" | awk '{print NF}') - 1))
echo "merged samples: ${n_samples}"

# level splits (header + level rows); species excludes SGB/strain (t__)
head -1 "${OUT}/metaphlan_merged_all.tsv" > "${OUT}/metaphlan_species.tsv"
awk -F'\t' '$1 ~ /\|s__/ && $1 !~ /\|t__/' "${OUT}/metaphlan_merged_all.tsv" >> "${OUT}/metaphlan_species.tsv"
head -1 "${OUT}/metaphlan_merged_all.tsv" > "${OUT}/metaphlan_genus.tsv"
awk -F'\t' '$1 ~ /\|g__/ && $1 !~ /\|s__/' "${OUT}/metaphlan_merged_all.tsv" >> "${OUT}/metaphlan_genus.tsv"

echo "=== classified summary ==="
python pipeline/04_taxonomy-metaphlan/metaphlan_classified_summary.py \
    "${OUT}/metaphlan_merged_all.tsv" "$ACC" "${OUT}/metaphlan_classified_summary.tsv"

echo "species rows: $(($(wc -l < "${OUT}/metaphlan_species.tsv") - 1))"
echo "genus rows:   $(($(wc -l < "${OUT}/metaphlan_genus.tsv") - 1))"
echo "mpa merge DONE $(date +'%F %H:%M:%S')"
