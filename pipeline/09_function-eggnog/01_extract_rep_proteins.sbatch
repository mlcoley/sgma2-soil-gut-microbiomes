#!/bin/bash
#SBATCH --job-name=reprot
#SBATCH --output=logs/reprot_%j.out
#SBATCH --error=logs/reprot_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Extracts the Prodigal protein sequence for every catalog representative by exact
# id, so eggNOG runs protein search (blastp) instead of translated nucleotide
# search (blastx), then splits the proteins into 250,000-sequence chunks for the
# eggNOG array.
# Output: genecat/catalog_rep_proteins.faa (one protein per catalog rep) and
#   genecat/eggnog_prot_chunks/catalog_rep_proteins.part_*.faa

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
echo "start $(date +'%F %H:%M:%S') on $(hostname)"

set +u
source ~/.bashrc
conda activate genecat
set -u

BASE="${SGMA2_BASE}"
WORK="${BASE}/assembly"
GC="${WORK}/genecat"
ORFS="${WORK}/orfs"
mkdir -p logs

REPFA="${GC}/catalog_rep_seq.fasta"
IDS="${GC}/rep_ids.txt"
POOL="${GC}/all_proteins.faa"
OUT="${GC}/catalog_rep_proteins.faa"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

CHUNKS="${GC}/eggnog_prot_chunks"

[ -f "$REPFA" ] || { echo "ERROR: catalog missing: $REPFA"; exit 1; }

# 1) rep ids (bare, first field) from the nucleotide catalog
echo "[1/5] extracting rep ids ... $(date +'%H:%M:%S')"
seqkit seq -n -i "$REPFA" > "$IDS"
nids=$(wc -l < "$IDS")
echo "    rep ids: $nids"

# 2) pool every sample's Prodigal proteins (ids already sample-prefixed + unique)
echo "[2/5] pooling orfs/*/*.faa ... $(date +'%H:%M:%S')"
cat "${ORFS}"/*/*.faa > "$POOL"
npool=$(grep -c '^>' "$POOL" || true)
echo "    pooled proteins: $npool"

# 3) pull the rep proteins by exact id
echo "[3/5] seqkit grep rep proteins ... $(date +'%H:%M:%S')"
seqkit grep -j "$THREADS" -f "$IDS" "$POOL" -o "$OUT"

# 4) verify count matches catalog
echo "[4/5] verifying count ... $(date +'%H:%M:%S')"
nout=$(grep -c '^>' "$OUT" || true)
echo "    rep proteins written: $nout   (expected $nids)"
[ "$nout" = "$nids" ] || { echo "ERROR: count mismatch ($nout != $nids) — some reps have no protein match"; exit 1; }

# 5) split into fixed-size chunks for eggNOG array
echo "[5/5] splitting into chunks ... $(date +'%H:%M:%S')"
seqkit split2 -s 250000 -O "$CHUNKS" -j "$THREADS" "$OUT"

# tidy the big intermediate pool. rep proteins + ids are kept
rm -f "$POOL"
echo "done $(date +'%F %H:%M:%S')  ->  $OUT"
