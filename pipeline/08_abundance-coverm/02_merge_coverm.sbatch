#!/bin/bash
#SBATCH --job-name=coverm_merge
#SBATCH --output=logs/coverm_merge_%j.out
#SBATCH --error=logs/coverm_merge_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=6:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Stacks the per-sample CoverM tables into one sparse long table (non-zero counts
# only): gene_id, sample_id, count, covered_fraction. Refuses to run until every
# manifest sample has a done marker.
# Output: assembly/genecat/gene_abundance_long.tsv.gz (no header),
#   coverm_samples.txt, coverm_merge_report.txt

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
set +u; source ~/.bashrc; conda activate genecat; set -u

BASE="${SGMA2_BASE}"
WORK="${BASE}/assembly"
GENECAT="${WORK}/genecat"
DONE="${WORK}/done_coverm"
MANIFEST="${WORK}/assembly_manifest.tsv"
OUT="${GENECAT}/gene_abundance_long.tsv.gz"
mkdir -p logs "$GENECAT"

nexp=$(wc -l < "$MANIFEST")
ndone=$(ls "$DONE"/*.done 2>/dev/null | wc -l)
echo "manifest samples: $nexp | done markers: $ndone"
[ "$ndone" -eq "$nexp" ] || { echo "ERROR: only $ndone/$nexp CoverM samples done — CoverM array not complete. Aborting."; exit 1; }

PIGZ="pigz -p ${SLURM_CPUS_PER_TASK:-8}"
command -v pigz >/dev/null 2>&1 || PIGZ="gzip"
echo "compressor: $PIGZ"

echo "Streaming merge -> $OUT   $(date +'%F %H:%M:%S')"
python3 pipeline/08_abundance-coverm/merge_coverm.py | $PIGZ > "${OUT}.partial"
mv "${OUT}.partial" "$OUT"

echo "=== output non-zero cells ==="
$PIGZ -dc "$OUT" | wc -l
echo "=== report head ==="
head -20 "${GENECAT}/coverm_merge_report.txt"
echo "Done $(date +'%F %H:%M:%S')"
