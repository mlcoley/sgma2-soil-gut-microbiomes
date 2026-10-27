#!/bin/bash
#SBATCH --job-name=amr
#SBATCH --output=logs/amr_%A_%a.out
#SBATCH --error=logs/amr_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# AMRFinderPlus screen of the catalog protein chunks (the same chunks used for
# eggNOG), protein mode with --plus and no --organism, so only acquired resistance
# genes plus stress, biocide and virulence elements are reported. Hits are keyed on
# the sample-prefixed catalog gene id. One array task per chunk: submit with
# --array=1-<number of chunk files>%<concurrency>.
# Output: assembly/genecat/amr/amrfinder/chunk_NNN.tsv

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
echo "Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} start $(date +'%F %H:%M:%S') on $(hostname)"

set +u; source ~/.bashrc; conda activate amr; set -u

BASE="${SGMA2_BASE}"
WORK="${BASE}/assembly"
CHUNKDIR="${WORK}/genecat/eggnog_prot_chunks"
DB="${WORK}/genecat/amr/amrfinder_db"
OUT="${WORK}/genecat/amr/amrfinder"
DONE="${WORK}/done_amr"
mkdir -p "$OUT" "$DONE" logs

# resolve the database written by 00_setup/04_amr_env_setup.sbatch
DBDIR="${DB}/latest"
[ -e "$DBDIR" ] || DBDIR=$(ls -d "${DB}"/*/ 2>/dev/null | sort | tail -1)
[ -n "$DBDIR" ] && [ -e "$DBDIR" ] || { echo "ERROR: no AMRFinderPlus DB under $DB — run 00_setup/04_amr_env_setup.sbatch first."; exit 1; }

[ -d "$CHUNKDIR" ] || { echo "ERROR: chunk dir missing: $CHUNKDIR"; exit 1; }
mapfile -t PARTS < <(ls "${CHUNKDIR}"/catalog_rep_proteins.part_*.faa 2>/dev/null | sort)
NPARTS="${#PARTS[@]}"
[ "$NPARTS" -gt 0 ] || { echo "ERROR: no part files in $CHUNKDIR"; exit 1; }

idx="$SLURM_ARRAY_TASK_ID"
[ "$idx" -le "$NPARTS" ] || { echo "Task $idx > $NPARTS parts — past end, exit."; exit 0; }
CHUNK="${PARTS[$((idx-1))]}"
tag=$(printf 'chunk_%03d' "$idx")
echo "Chunk $idx/$NPARTS -> $(basename "$CHUNK")  (tag $tag)  DB $DBDIR"

[ -f "${DONE}/${tag}.done" ] && { echo "$tag already done — exit."; exit 0; }
[ -f "$CHUNK" ] || { echo "ERROR: chunk file missing: $CHUNK"; exit 1; }

TMP="${WORK}/tmp_amr/${tag}.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
rm -rf "$TMP"; mkdir -p "$TMP"
export TMPDIR="$TMP"
trap 'rm -rf "$TMP"' EXIT

THREADS="${SLURM_CPUS_PER_TASK:-16}"

echo "amrfinder -p --plus (protein, no --organism) ... $(date +'%H:%M:%S')"
amrfinder \
    -p "$CHUNK" \
    --plus \
    -d "$DBDIR" \
    --threads "$THREADS" \
    -o "${TMP}/${tag}.tsv"

[ -f "${TMP}/${tag}.tsv" ] || { echo "ERROR: no output produced for $tag"; exit 1; }
nhit=$(grep -vc '^Protein identifier' "${TMP}/${tag}.tsv" || true)
echo "=== $tag: ${nhit} ARG hits ==="

# atomic publish
mv "${TMP}/${tag}.tsv" "${OUT}/${tag}.tsv.partial"
mv "${OUT}/${tag}.tsv.partial" "${OUT}/${tag}.tsv"

touch "${DONE}/${tag}.done"
echo "Done $tag at $(date +'%F %H:%M:%S')"
