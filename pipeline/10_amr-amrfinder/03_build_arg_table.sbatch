#!/bin/bash
#SBATCH --job-name=arg_table
#SBATCH --output=logs/arg_table_%j.out
#SBATCH --error=logs/arg_table_%j.err
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
# Aggregates gene abundance (gene_abundance_long.tsv.gz) to ARG gene symbol using
# arg_hits_long.tsv.gz. Set MIN_COVERED_FRACTION (for example 0.5) to keep only gene
# x sample cells at or above that breadth; outputs then carry a _cf50 suffix.
# Output: genecat/amr/arg_abundance_long.tsv.gz, arg_abundance_wide.tsv.gz,
#   arg_annotation.tsv, arg_table_report.txt

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
set +u; source ~/.bashrc; conda activate genecat; set -u

BASE="${SGMA2_BASE}"
GEN="${BASE}/assembly/genecat"
mkdir -p logs

for f in "${GEN}/amr/arg_hits_long.tsv.gz" "${GEN}/gene_abundance_long.tsv.gz" "${GEN}/coverm_samples.txt"; do
  [ -s "$f" ] || { echo "ERROR: missing input $f"; exit 1; }
done

python pipeline/10_amr-amrfinder/build_arg_table.py

echo "==== report ===="
cat "${GEN}/amr/arg_table_report.txt"
