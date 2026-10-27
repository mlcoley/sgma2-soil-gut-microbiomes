#!/bin/bash
#SBATCH --job-name=breadth_rebuild
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/breadth_rebuild_%j.out
#SBATCH --error=logs/breadth_rebuild_%j.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Rebuilds the KO and ARG abundance tables keeping only gene x sample cells with
# CoverM covered_fraction >= 0.5, and writes the ARG gene slice with covered_fraction.
# Outputs carry a _cf50 suffix and sit alongside the unfiltered tables, giving a
# breadth-based sensitivity analysis of the presence-based results.
# Output: genecat/eggnog/ko_abundance_wide_cf50.tsv.gz,
#   genecat/amr/arg_abundance_wide_cf50.tsv.gz, genecat/amr/arg_gene_abundance.tsv.gz

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
echo "breadth rebuild job ${SLURM_JOB_ID} start $(date +'%F %H:%M:%S') on $(hostname)"

set +u
source ~/.bashrc
conda activate genecat
set -u

BASE="${SGMA2_BASE}"
WORK="${BASE}/assembly"
GC="${WORK}/genecat"
export MIN_COVERED_FRACTION=0.5
mkdir -p logs

for f in "$GC/gene_abundance_long.tsv.gz" "$GC/coverm_samples.txt" \
         "$GC/eggnog/catalog_ko_ec.tsv.gz" "$GC/amr/arg_hits_long.tsv.gz"; do
  [ -s "$f" ] || { echo "FATAL: missing input $f" >&2; exit 1; }
done

echo "=== 1. KO wide table at covered_fraction >= ${MIN_COVERED_FRACTION} ==="
python pipeline/09_function-eggnog/build_ko_table.py
echo "--- KO report:"; cat "$GC/eggnog/ko_table_report_cf50.txt" || true

echo
echo "=== 2. ARG wide table at covered_fraction >= ${MIN_COVERED_FRACTION} ==="
python pipeline/10_amr-amrfinder/build_arg_table.py
echo "--- ARG report:"; cat "$GC/amr/arg_table_report_cf50.txt" || true

echo
echo "=== 3. ARG gene (allele) slice, now KEEPING covered_fraction ==="
bash pipeline/11_feature-tables/extract_arg_gene_abundance.sbatch

echo
echo "=== DONE ==="
ls -la "$GC/eggnog/ko_abundance_wide_cf50.tsv.gz" \
       "$GC/amr/arg_abundance_wide_cf50.tsv.gz" \
       "$GC/amr/arg_gene_abundance.tsv.gz"
