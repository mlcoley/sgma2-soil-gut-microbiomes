#!/bin/bash
#SBATCH --job-name=eggnog_merge
#SBATCH --output=logs/eggnog_merge_%j.out
#SBATCH --error=logs/eggnog_merge_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Concatenates the per-chunk eggNOG annotation files (one header, no '##' lines)
# and extracts a slim table of gene_id, KEGG_ko and EC. Does not run until every
# chunk has a done marker. KEGG_ko and EC can be comma-separated lists.
# Output: genecat/eggnog/catalog_eggnog_annotations.tsv.gz, catalog_ko_ec.tsv.gz,
#   eggnog_merge_report.txt

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
set +u; source ~/.bashrc; conda activate genecat; set -u

BASE="${SGMA2_BASE}"
WORK="${BASE}/assembly"
EGG="${WORK}/genecat/eggnog"
DONE="${WORK}/done_eggnog"
CHUNKDIR="${WORK}/genecat/eggnog_prot_chunks"
FULL="${EGG}/catalog_eggnog_annotations.tsv.gz"
SLIM="${EGG}/catalog_ko_ec.tsv.gz"
REPORT="${EGG}/eggnog_merge_report.txt"
mkdir -p logs "$EGG"

# ---- gate on the part count: one done-marker per chunk ----
NPARTS=$(ls "${CHUNKDIR}"/catalog_rep_proteins.part_*.faa 2>/dev/null | wc -l)
ndone=$(ls "$DONE"/chunk_*.done 2>/dev/null | wc -l)
echo "chunks expected: $NPARTS | done markers: $ndone"
[ "$NPARTS" -gt 0 ] || { echo "ERROR: no chunk part files found in $CHUNKDIR"; exit 1; }
[ "$ndone" -eq "$NPARTS" ] || { echo "ERROR: only $ndone/$NPARTS eggNOG chunks done — eggNOG array not complete. Aborting."; exit 1; }

mapfile -t ANN < <(ls "${EGG}"/chunk_*.emapper.annotations 2>/dev/null | sort)
nann="${#ANN[@]}"
echo "annotation files: $nann"
[ "$nann" -eq "$NPARTS" ] || { echo "ERROR: $nann annotation files but $NPARTS chunks — missing outputs. Aborting."; exit 1; }

PIGZ="pigz -p ${SLURM_CPUS_PER_TASK:-8}"
command -v pigz >/dev/null 2>&1 || PIGZ="gzip"
echo "compressor: $PIGZ"

# ---- full table: one '#query' header (from the first chunk), all data rows, no '##' ----
echo "Merging full annotations -> $FULL   $(date +'%F %H:%M:%S')"
{
  grep -m1 '^#query' "${ANN[0]}"
  for f in "${ANN[@]}"; do grep -v '^#' "$f"; done
} | $PIGZ > "${FULL}.partial"
mv "${FULL}.partial" "$FULL"

# ---- slim gene_id -> KO/EC table (cols 1,12,11) ----
echo "Building slim KO/EC table -> $SLIM   $(date +'%F %H:%M:%S')"
{
  printf 'gene_id\tKEGG_ko\tEC\n'
  for f in "${ANN[@]}"; do awk -F'\t' '!/^#/{print $1"\t"$12"\t"$11}' "$f"; done
} | $PIGZ > "${SLIM}.partial"
mv "${SLIM}.partial" "$SLIM"

# ---- report ----
echo "Writing report -> $REPORT"
{
  echo "eggNOG merge report  $(date +'%F %H:%M:%S')"
  echo "chunks merged: $nann / $NPARTS"
  tot=$($PIGZ -dc "$FULL" | grep -vc '^#' || true)
  ko=$($PIGZ -dc "$SLIM" | awk -F'\t' 'NR>1 && $2!="-"' | wc -l)
  ec=$($PIGZ -dc "$SLIM" | awk -F'\t' 'NR>1 && $3!="-"' | wc -l)
  uq=$($PIGZ -dc "$SLIM" | awk -F'\t' 'NR>1{print $1}' | sort -u | wc -l)
  echo "annotated genes (data rows): $tot"
  echo "  with KEGG_ko (!= '-'):     $ko"
  echo "  with EC (!= '-'):          $ec"
  echo "  unique gene_ids:           $uq"
  echo "catalog size (for context):  $(wc -l < "${WORK}/genecat/catalog_gene_lengths.tsv") reps"
} | tee "$REPORT"

echo "Done $(date +'%F %H:%M:%S')"
