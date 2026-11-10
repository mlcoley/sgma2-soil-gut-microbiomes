#!/bin/bash
#SBATCH --job-name=mmseqs_cat
#SBATCH --output=logs/mmseqs_catalog_%j.out
#SBATCH --error=logs/mmseqs_catalog_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=12:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Pools all per-sample gene nucleotide sequences (.ffn) and builds a non-redundant
# gene catalog with mmseqs easy-linclust: 95% nucleotide identity, 90% coverage of
# the shorter sequence (--cov-mode 1). Version: mmseqs2 18.8cc5c.
# Output: assembly/genecat/ catalog_rep_seq.fasta, catalog_cluster.tsv,
#   catalog_all_seqs.fasta, catalog_stats.txt, catalog_gene_lengths.tsv

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
echo "mmseqs catalog job ${SLURM_JOB_ID} start $(date +'%F %H:%M:%S') on $(hostname)"

set +u
source ~/.bashrc
conda activate genecat
set -u

BASE="${SGMA2_BASE}"
WORK="${BASE}/assembly"
ORFS="${WORK}/orfs"
OUT="${WORK}/genecat"
TMP="${WORK}/tmp_mmseqs.${SLURM_JOB_ID}"
mkdir -p "$OUT" logs
rm -rf "$TMP"; mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

THREADS="${SLURM_CPUS_PER_TASK:-32}"

# ---- 1. pool every sample's nucleotide genes into one fasta ----
POOL="${OUT}/all_genes.ffn"
echo "Pooling .ffn from $(ls -d ${ORFS}/*/ | wc -l) samples ... $(date +'%H:%M:%S')"
: > "$POOL"
nsamp=0
for d in "${ORFS}"/*/; do
    s=$(basename "$d")
    ffn="${d}/${s}.ffn"
    [ -f "$ffn" ] || { echo "ERROR: missing $ffn"; exit 1; }
    cat "$ffn" >> "$POOL"
    nsamp=$((nsamp+1))
done
n_in=$(grep -c '^>' "$POOL")
echo "Pooled ${nsamp} samples -> ${n_in} genes into ${POOL}"

# ---- 2. linclust non-redundant catalog ----
echo "mmseqs easy-linclust (95% id / 90% cov-mode 1) ... $(date +'%H:%M:%S')"
mmseqs easy-linclust \
    "$POOL" \
    "${OUT}/catalog" \
    "$TMP" \
    --min-seq-id 0.95 \
    -c 0.9 \
    --cov-mode 1 \
    --kmer-per-seq 80 \
    --threads "$THREADS" \
    --split-memory-limit 200G

# ---- 3. stats ----
n_rep=$(grep -c '^>' "${OUT}/catalog_rep_seq.fasta")
{
    echo "SGMA2 gene catalog — $(date +'%F %H:%M:%S')"
    echo "samples pooled      : ${nsamp}"
    echo "input genes         : ${n_in}"
    echo "catalog reps (NR)   : ${n_rep}"
    awk -v i="$n_in" -v r="$n_rep" 'BEGIN{printf "redundancy collapse : %.2f%% (reps / input)\n", 100*r/i}'
    echo "params              : --min-seq-id 0.95 -c 0.9 --cov-mode 1 (linclust, nucleotide)"
} | tee "${OUT}/catalog_stats.txt"

# ---- 4. gene length table (gene_id, length), used for length normalisation ----
seqkit fx2tab -nli "${OUT}/catalog_rep_seq.fasta" | cut -f1,2 > "${OUT}/catalog_gene_lengths.tsv"

echo "Done catalog: ${n_rep} representatives from ${n_in} genes at $(date +'%F %H:%M:%S')"
