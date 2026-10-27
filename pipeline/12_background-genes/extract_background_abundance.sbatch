#!/bin/bash
#SBATCH --job-name=bg_abundance
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/bg_abundance.%j.out
#SBATCH --error=logs/bg_abundance.%j.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Extracts per-sample abundance rows for the selected background genes
# (GENES, default background_gene_ids.txt under $SGMA2_BASE/assembly). Set TAG to write
# a supplementary extraction without re-streaming genes already pulled:
#   TAG=extra GENES=<ids file> sbatch pipeline/12_background-genes/extract_background_abundance.sbatch
# Output: assembly/background/background_abundance[_TAG].tsv.gz

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail

set +u
source ~/.bashrc
conda activate genecat
set -u

BASE=${SGMA2_BASE}/assembly
GC=$BASE/genecat
OUTDIR=$BASE/background

TAG=${TAG:-}
GENES=${GENES:-background_gene_ids.txt}
OUT=$OUTDIR/background_abundance${TAG:+_$TAG}.tsv.gz
DONE=$OUTDIR/done/job2${TAG:+_$TAG}.done

mkdir -p "$OUTDIR" "$OUTDIR/done" "$BASE/logs"

if [[ -f "$DONE" ]]; then
  echo "job2 already complete ($DONE) — nothing to do."
  exit 0
fi

TMP=$(mktemp -d "$BASE/tmp_bg2.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

echo "host  : $(hostname)"
echo "start : $(date -Is)"

for f in "$GC/gene_abundance_long.tsv.gz" \
         "$BASE/T23c_cohort_compartments.tsv" \
         "$BASE/$GENES"; do
  [[ -s "$f" ]] || { echo "MISSING OR EMPTY: $f" >&2; exit 1; }
done

python pipeline/12_background-genes/extract_background_abundance.py \
  --abundance "$GC/gene_abundance_long.tsv.gz" \
  --genes     "$BASE/$GENES" \
  --cohort    "$BASE/T23c_cohort_compartments.tsv" \
  --out       "$TMP/out.tsv.gz"

mv "$TMP/out.tsv.gz" "$OUT"
touch "$DONE"

echo "genes : $GENES"
echo "out   : $OUT"
echo "rows  : $(zcat "$OUT" | wc -l)"
echo "size  : $(du -h "$OUT" | cut -f1)"
echo "end   : $(date -Is)"
