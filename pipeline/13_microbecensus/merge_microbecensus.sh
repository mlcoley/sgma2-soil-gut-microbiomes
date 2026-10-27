#!/bin/bash
# Assemble per-sample MicrobeCensus outputs -> feature_tables/microbecensus_ge.tsv
# Run after the array finishes; then pull the TSV local into feature_tables/.
# Usage:  bash merge_microbecensus.sh  <per_sample_dir>  <out.tsv>
set -euo pipefail
DIR=${1:-per_sample}
OUT=${2:-microbecensus_ge.tsv}
printf 'sample_id\taverage_genome_size\tgenome_equivalents\n' > "$OUT"
for f in "$DIR"/*.mc.txt; do
  s=$(basename "$f" .mc.txt)
  ags=$(awk -F'\t' '/^average_genome_size:/{print $2}' "$f")
  ge=$(awk -F'\t'  '/^genome_equivalents:/{print $2}'  "$f")
  printf '%s\t%s\t%s\n' "$s" "$ags" "$ge" >> "$OUT"
done
echo "wrote $OUT ($(($(wc -l < "$OUT") - 1)) samples)"
