#!/bin/bash
#SBATCH --job-name=eggnog
#SBATCH --output=logs/eggnog_%A_%a.out
#SBATCH --error=logs/eggnog_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# eggNOG-mapper 2.1.13 annotation (diamond 2.2.2, eggNOG 5.0) of the catalog protein
# chunks, giving KO and EC assignments keyed on catalog gene id. One array task per
# chunk: submit with --array=1-<number of chunk files>%<concurrency>.
# The shared database (eggnog.db, eggnog_proteins.dmnd) is read from
# $SGMA2_DB_DIR/eggnog and is read-only. --dmnd_iterate no is required for the
# format-v2 diamond index; --dbmem loads eggnog.db into memory (about 44 GB, so
# --mem=96G).
# Output: assembly/genecat/eggnog/chunk_NNN.emapper.{annotations,hits,seed_orthologs}

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
echo "Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} start $(date +'%F %H:%M:%S') on $(hostname)"

set +u
source ~/.bashrc
conda activate genecat
set -u

BASE="${SGMA2_BASE}"
WORK="${BASE}/assembly"
EGGDB="${SGMA2_DB_DIR}/eggnog"
CHUNKDIR="${WORK}/genecat/eggnog_prot_chunks"
OUT="${WORK}/genecat/eggnog"
DONE="${WORK}/done_eggnog"
mkdir -p "$OUT" "$DONE" logs

[ -d "$EGGDB" ]                        || { echo "ERROR: shared eggNOG DB dir missing: $EGGDB"; exit 1; }
[ -f "${EGGDB}/eggnog.db" ]            || { echo "ERROR: eggnog.db missing in $EGGDB"; exit 1; }
[ -f "${EGGDB}/eggnog_proteins.dmnd" ] || { echo "ERROR: eggnog_proteins.dmnd missing in $EGGDB"; exit 1; }

# ---- this task's chunk: the Nth part file in sorted order (padding-agnostic) ----
[ -d "$CHUNKDIR" ] || { echo "ERROR: chunk dir missing: $CHUNKDIR — run the seqkit split2 step first (see header)."; exit 1; }
mapfile -t PARTS < <(ls "${CHUNKDIR}"/catalog_rep_proteins.part_*.faa 2>/dev/null | sort)
NPARTS="${#PARTS[@]}"
[ "$NPARTS" -gt 0 ] || { echo "ERROR: no part files in $CHUNKDIR — split step not run?"; exit 1; }

idx="$SLURM_ARRAY_TASK_ID"
[ "$idx" -le "$NPARTS" ] || { echo "Task $idx > $NPARTS parts — past end, exit."; exit 0; }
CHUNK="${PARTS[$((idx-1))]}"                 # 0-based array, 1-based task id
tag=$(printf 'chunk_%03d' "$idx")           # stable output tag regardless of seqkit's padding
echo "Chunk $idx/$NPARTS -> $(basename "$CHUNK")  (tag $tag)"

[ -f "${DONE}/${tag}.done" ] && { echo "$tag already done — exit."; exit 0; }
[ -f "$CHUNK" ] || { echo "ERROR: chunk file missing: $CHUNK"; exit 1; }

# job-private workspace keyed by the SLURM job id -> two jobs on the same chunk
# can never collide; emapper's temp + partial outputs live here, discarded on exit.
TMP="${WORK}/tmp_eggnog/${tag}.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
rm -rf "$TMP"; mkdir -p "${TMP}/out"
export TMPDIR="$TMP"
trap 'rm -rf "$TMP"' EXIT

THREADS="${SLURM_CPUS_PER_TASK:-16}"

echo "emapper (diamond, --itype proteins, --dmnd_iterate no) ... $(date +'%H:%M:%S')"
emapper.py \
    -m diamond --itype proteins \
    -i "$CHUNK" \
    --data_dir "$EGGDB" \
    --dmnd_db "${EGGDB}/eggnog_proteins.dmnd" \
    --dmnd_iterate no \
    --dbmem \
    --output "$tag" --output_dir "${TMP}/out" \
    --temp_dir "$TMP" \
    --cpu "$THREADS" \
    --override

ANN="${TMP}/out/${tag}.emapper.annotations"
[ -f "$ANN" ] || { echo "ERROR: no annotations produced for $tag"; exit 1; }
nann=$(grep -vc '^#' "$ANN" || true)
echo "=== $tag: ${nann} annotated genes ==="

# atomic publish: only a fully-successful run reaches ${OUT}
for ext in emapper.annotations emapper.hits emapper.seed_orthologs; do
    f="${TMP}/out/${tag}.${ext}"
    [ -f "$f" ] && { mv "$f" "${OUT}/${tag}.${ext}.partial"; mv "${OUT}/${tag}.${ext}.partial" "${OUT}/${tag}.${ext}"; }
done

touch "${DONE}/${tag}.done"
echo "Done $tag at $(date +'%F %H:%M:%S')"
