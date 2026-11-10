#!/bin/bash
#SBATCH --job-name=prodigal
#SBATCH --output=logs/prodigal_%A_%a.out
#SBATCH --error=logs/prodigal_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Per-assembly gene calling with Prodigal (-p meta); partial genes are kept.
# Gene and contig ids are prefixed with the sample id (k141_123_1 -> P###_k141_123_1)
# because every MEGAHIT assembly reuses the same contig names and the next step pools
# all genes into one catalog.
# One array task per row of assembly_manifest.tsv; set --array=1-N.
# Output: assembly/orfs/<sample>/<sample>.{faa,ffn,gff}

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
echo "Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} start $(date +'%F %H:%M:%S') on $(hostname)"

set +u
source ~/.bashrc
conda activate genecat
set -u

BASE="${SGMA2_BASE}"
WORK="${BASE}/assembly"
ASM="${WORK}/assemblies"
OUT="${WORK}/orfs"
DONE="${WORK}/done_orf"
MANIFEST="${WORK}/assembly_manifest.tsv"
mkdir -p "$OUT" "$DONE" logs

# ---- this task's (sample, compartment) from the manifest ----
[ -f "$MANIFEST" ] || { echo "ERROR: manifest missing: $MANIFEST"; exit 1; }
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
[ -z "$line" ] && { echo "No row ${SLURM_ARRAY_TASK_ID} — past end of manifest, exit."; exit 0; }
s=$(echo "$line"    | cut -f1)
comp=$(echo "$line" | cut -f2)
echo "Sample: $s   Compartment: $comp"
[ -f "${DONE}/${s}.done" ] && { echo "$s already done — exit."; exit 0; }

CONTIGS_IN="${ASM}/${s}/${s}.contigs.fa"
[ -f "$CONTIGS_IN" ] || { echo "ERROR: contigs missing: $CONTIGS_IN (assembly step incomplete for $s)"; exit 1; }

# job-private workspace keyed by the SLURM job id -> two concurrent jobs on the
# same sample can never share a path.
TMP="${WORK}/tmp_orf/${s}.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
rm -rf "$TMP"; mkdir -p "${TMP}/out"
trap 'rm -rf "$TMP"' EXIT

echo "Prodigal (-p meta) ... $(date +'%H:%M:%S')"
prodigal -i "$CONTIGS_IN" -p meta -f gff -q \
         -a "${TMP}/raw.faa" -d "${TMP}/raw.ffn" -o "${TMP}/raw.gff"

# ---- prefix sample id onto every gene/contig id (see note above) ----
sed "s/^>/>${s}_/" "${TMP}/raw.faa" > "${TMP}/out/${s}.faa"
sed "s/^>/>${s}_/" "${TMP}/raw.ffn" > "${TMP}/out/${s}.ffn"
awk -v s="$s" 'BEGIN{FS=OFS="\t"} /^#/{print;next} NF==0{print;next} {$1=s"_"$1; print}' \
    "${TMP}/raw.gff" > "${TMP}/out/${s}.gff"

ngenes=$(grep -c '^>' "${TMP}/out/${s}.faa" || true)
echo "=== $s ($comp): ${ngenes} genes called ==="

# atomic publish: only a fully-successful run reaches ${OUT}/${s}
rm -rf "${OUT}/${s}.partial" "${OUT}/${s}"
mv "${TMP}/out" "${OUT}/${s}.partial"
mv "${OUT}/${s}.partial" "${OUT}/${s}"

touch "${DONE}/${s}.done"
echo "Done $s (${ngenes} genes) at $(date +'%F %H:%M:%S')"
