#!/bin/bash
#SBATCH --job-name=bg_gene_summary
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/bg_gene_summary.%j.out
#SBATCH --error=logs/bg_gene_summary.%j.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Reduces gene_abundance_long.tsv.gz to one row per detected catalog gene:
# per-compartment prevalence under two presence rules (count >= 1 and covered
# fraction >= 0.5), abundance, gene length, eggNOG orthologous group and KO. Genes
# called by AMRFinderPlus and genes without an eggNOG annotation are excluded.
# Needs $SGMA2_BASE/assembly/{T23c_cohort_compartments.tsv,T23d_amrfinder_gene_ids.txt},
# written by analysis/util_background_target_profile.R.
# Also writes arg_gene_lengths.tsv (lengths of the AMRFinderPlus genes).
# Output: assembly/background/background_gene_summary.tsv.gz

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail

set +u
source ~/.bashrc
conda activate genecat
set -u

BASE=${SGMA2_BASE}/assembly
GC=$BASE/genecat
OUTDIR=$BASE/background
DONE=$OUTDIR/done/job1.done

mkdir -p "$OUTDIR" "$OUTDIR/done" "$BASE/logs"

if [[ -f "$DONE" ]]; then
  echo "job1 already complete ($DONE) — nothing to do."
  exit 0
fi

TMP=$(mktemp -d "$BASE/tmp_bg.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

echo "host        : $(hostname)"
echo "start       : $(date -Is)"
echo "python      : $(python --version 2>&1)"
echo "abundance   : $GC/gene_abundance_long.tsv.gz"
echo "lengths     : $GC/catalog_gene_lengths.tsv"
echo "eggnog      : $GC/eggnog/catalog_eggnog_annotations.tsv.gz"

for f in "$GC/gene_abundance_long.tsv.gz" \
         "$GC/catalog_gene_lengths.tsv" \
         "$GC/eggnog/catalog_eggnog_annotations.tsv.gz" \
         "$BASE/T23c_cohort_compartments.tsv" \
         "$BASE/T23d_amrfinder_gene_ids.txt"; do
  [[ -s "$f" ]] || { echo "MISSING OR EMPTY: $f" >&2; exit 1; }
done

python pipeline/12_background-genes/background_gene_summary.py \
  --abundance "$GC/gene_abundance_long.tsv.gz" \
  --lengths   "$GC/catalog_gene_lengths.tsv" \
  --eggnog    "$GC/eggnog/catalog_eggnog_annotations.tsv.gz" \
  --cohort    "$BASE/T23c_cohort_compartments.tsv" \
  --exclude   "$BASE/T23d_amrfinder_gene_ids.txt" \
  --out       "$TMP/background_gene_summary.tsv.gz"

mv "$TMP/background_gene_summary.tsv.gz" "$OUTDIR/background_gene_summary.tsv.gz"
grep -F -w -f "$BASE/T23d_amrfinder_gene_ids.txt" "$GC/catalog_gene_lengths.tsv" > "$BASE/arg_gene_lengths.tsv"
touch "$DONE"

echo "rows  : $(zcat "$OUTDIR/background_gene_summary.tsv.gz" | wc -l)"
echo "size  : $(du -h "$OUTDIR/background_gene_summary.tsv.gz" | cut -f1)"
echo "end   : $(date -Is)"
