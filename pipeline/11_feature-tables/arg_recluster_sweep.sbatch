#!/bin/bash
#SBATCH --job-name=arg_reclust
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/arg_reclust_%j.out
#SBATCH --error=logs/arg_reclust_%j.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Re-clusters the ARG-bearing catalog genes (nucleotide) at 95, 90, 85 and 80%
# identity with mmseqs easy-linclust (-c 0.9 --cov-mode 1, as in the catalog build)
# for the allele-threshold sensitivity analysis. The sequences are already 95%
# representatives, so the sweep can only merge further; 0.95 is a control.
# Output: assembly/genecat/arg_reclust/arg_id{95,90,85,80}_cluster.tsv

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
echo "arg reclust job ${SLURM_JOB_ID} start $(date +'%F %H:%M:%S') on $(hostname)"

set +u
source ~/.bashrc
conda activate genecat
set -u

BASE="${SGMA2_BASE}"
WORK="${BASE}/assembly"
GC="${WORK}/genecat"
OUT="${GC}/arg_reclust"
TMP="${WORK}/tmp_arg_reclust.${SLURM_JOB_ID}"
mkdir -p "$OUT" "$TMP" logs
trap 'rm -rf "$TMP"' EXIT

for f in "$GC/amr/arg_hits_long.tsv.gz" "$GC/catalog_rep_seq.fasta"; do
  [ -s "$f" ] || { echo "FATAL: missing input $f" >&2; exit 1; }
done
echo "inputs found under $GC"

echo "=== 1. ARG gene IDs from the AMR merge ==="
zcat "$GC/amr/arg_hits_long.tsv.gz" | awk -F'\t' 'NR>1 && $1!="" {print $1}' | sort -u > "$OUT/arg_gene_ids.txt"
N_IDS=$(wc -l < "$OUT/arg_gene_ids.txt")
echo "ARG-bearing catalog genes: $N_IDS"

echo "=== 2. pull those sequences out of the catalog ==="
seqkit grep -f "$OUT/arg_gene_ids.txt" "$GC/catalog_rep_seq.fasta" > "$OUT/arg_rep_seqs.fasta"
N_SEQ=$(grep -c '^>' "$OUT/arg_rep_seqs.fasta")
echo "sequences extracted: $N_SEQ"
if [ "$N_SEQ" -ne "$N_IDS" ]; then
  echo "FATAL: extracted $N_SEQ sequences for $N_IDS ids — aborting rather than sweeping an incomplete set" >&2
  exit 1
fi
seqkit stats "$OUT/arg_rep_seqs.fasta"

echo "=== 3. re-cluster at each identity threshold ==="
for ID in 0.95 0.90 0.85 0.80; do
  TAG=$(echo "$ID" | sed 's/0\.//')          # 95, 90, 85, 80
  echo "--- --min-seq-id $ID"
  rm -rf "${TMP:?}/mm_$TAG"; mkdir -p "$TMP/mm_$TAG"
  mmseqs easy-linclust "$OUT/arg_rep_seqs.fasta" "$OUT/arg_id${TAG}" "$TMP/mm_$TAG" \
      --min-seq-id "$ID" -c 0.9 --cov-mode 1 --threads "${SLURM_CPUS_PER_TASK:-8}" >/dev/null
  # easy-linclust writes <prefix>_cluster.tsv as  representative <TAB> member
  NCL=$(cut -f1 "$OUT/arg_id${TAG}_cluster.tsv" | sort -u | wc -l)
  NMEM=$(wc -l < "$OUT/arg_id${TAG}_cluster.tsv")
  echo "    id=$ID -> $NCL clusters covering $NMEM genes"
  rm -f "$OUT/arg_id${TAG}_all_seqs.fasta"   # large and not needed downstream
done

echo "=== DONE ==="
ls -la "$OUT"/arg_id*_cluster.tsv
