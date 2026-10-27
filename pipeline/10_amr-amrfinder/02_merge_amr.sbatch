#!/bin/bash
#SBATCH --job-name=amr_merge
#SBATCH --output=logs/amr_merge_%j.out
#SBATCH --error=logs/amr_merge_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Concatenates the per-chunk AMRFinderPlus tables into a native table and a
# harmonised hit table (arg_hits_long.tsv.gz) keyed on gene_id with a tool column.
# Refuses to run until every chunk has a done marker.
# Output: genecat/amr/amrfinder/catalog_amrfinder.tsv.gz, genecat/amr/arg_hits_long.tsv.gz,
#   genecat/amr/amr_merge_report.txt

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
set +u; source ~/.bashrc; conda activate amr; set -u

BASE="${SGMA2_BASE}"
WORK="${BASE}/assembly"
CHUNKDIR="${WORK}/genecat/eggnog_prot_chunks"
OUT="${WORK}/genecat/amr/amrfinder"
DONE="${WORK}/done_amr"
mkdir -p logs

NPARTS=$(ls "${CHUNKDIR}"/catalog_rep_proteins.part_*.faa 2>/dev/null | wc -l)
ndone=$(ls "$DONE"/chunk_*.done 2>/dev/null | wc -l)
nout=$(ls "$OUT"/chunk_*.tsv 2>/dev/null | wc -l)
echo "chunks expected: $NPARTS | done: $ndone | output tsvs: $nout"
[ "$NPARTS" -gt 0 ] || { echo "ERROR: no chunk parts found in $CHUNKDIR"; exit 1; }
[ "$ndone" -eq "$NPARTS" ] || { echo "ERROR: only $ndone/$NPARTS AMR chunks done — array not complete. Aborting."; exit 1; }
[ "$nout" -eq "$NPARTS" ] || { echo "ERROR: $nout/$NPARTS output tsvs present — missing chunk output. Aborting."; exit 1; }

python pipeline/10_amr-amrfinder/merge_amr.py

echo "==== report ===="
cat "${WORK}/genecat/amr/amr_merge_report.txt"
