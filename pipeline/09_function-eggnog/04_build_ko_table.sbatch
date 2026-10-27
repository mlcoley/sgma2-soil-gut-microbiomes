#!/bin/bash
#SBATCH --job-name=ko_table
#SBATCH --output=logs/ko_table_%j.out
#SBATCH --error=logs/ko_table_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Aggregates gene abundance (gene_abundance_long.tsv.gz) to KEGG Orthology using
# catalog_ko_ec.tsv.gz. Set MIN_COVERED_FRACTION (for example 0.5) to keep only gene
# x sample cells at or above that breadth; outputs then carry a _cf50 suffix.
# Output: genecat/eggnog/ko_abundance_long.tsv.gz, ko_abundance_wide.tsv.gz,
#   ko_table_report.txt

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
set +u; source ~/.bashrc; conda activate genecat; set -u

BASE="${SGMA2_BASE}"
GEN="${BASE}/assembly/genecat"
mkdir -p logs

for f in "${GEN}/eggnog/catalog_ko_ec.tsv.gz" "${GEN}/gene_abundance_long.tsv.gz" "${GEN}/coverm_samples.txt"; do
  [ -s "$f" ] || { echo "ERROR: missing input $f"; exit 1; }
done

python pipeline/09_function-eggnog/build_ko_table.py

echo "==== report ===="
cat "${GEN}/eggnog/ko_table_report.txt"
