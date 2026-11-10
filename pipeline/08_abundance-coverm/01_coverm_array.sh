#!/bin/bash
#SBATCH --job-name=coverm
#SBATCH --output=logs/coverm_%A_%a.out
#SBATCH --error=logs/coverm_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=8:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Per-sample read mapping to the gene catalog with CoverM (minimap2-sr). Produces raw
# read count and covered fraction for every catalog gene (--min-covered-fraction 0),
# so presence thresholds are applied downstream. Counts are not normalised here.
# One array task per row of assembly_manifest.tsv; submit with --array=1-N%8.
# SAMPLE=<id> runs a single sample by name. The control libraries are mapped to the
# same catalog with MANIFEST=<control_manifest.tsv> DONE=<dir> overrides.
# Version: coverm 0.7.0 (minimap2 2.31, samtools 1.23.1).
# Output: assembly/coverm/<sample>.tsv

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
echo "Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} start $(date +'%F %H:%M:%S') on $(hostname)"

set +u
source ~/.bashrc
conda activate genecat
set -u

BASE="${SGMA2_BASE}"
WORK="${BASE}/assembly"
READS="${BASE}/host_filter/filtered_fastq"
CATALOG="${WORK}/genecat/catalog_rep_seq.fasta"
OUT="${WORK}/coverm"
DONE="${DONE:-${WORK}/done_coverm}"
MANIFEST="${MANIFEST:-${WORK}/assembly_manifest.tsv}"
mkdir -p "$OUT" "$DONE" logs

[ -f "$CATALOG" ] || { echo "ERROR: catalog missing: $CATALOG (catalog step incomplete)"; exit 1; }

# ---- this task's sample: SAMPLE= env override (test) else manifest row ----
if [ -n "${SAMPLE:-}" ]; then
    s="$SAMPLE"
    comp=$(awk -F'\t' -v s="$s" '$1==s{print $2}' "$MANIFEST")
    echo "SAMPLE override = $s (test mode)"
else
    [ -f "$MANIFEST" ] || { echo "ERROR: manifest missing: $MANIFEST"; exit 1; }
    line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
    [ -z "$line" ] && { echo "No row ${SLURM_ARRAY_TASK_ID} — past end of manifest, exit."; exit 0; }
    s=$(echo "$line"    | cut -f1)
    comp=$(echo "$line" | cut -f2)
fi
echo "Sample: $s   Compartment: ${comp:-unknown}"
[ -f "${DONE}/${s}.done" ] && { echo "$s already done — exit."; exit 0; }

R1="${READS}/${s}_1.filt.fq.gz"
R2="${READS}/${s}_2.filt.fq.gz"
[ -f "$R1" ] || { echo "ERROR: R1 missing: $R1"; exit 1; }
[ -f "$R2" ] || { echo "ERROR: R2 missing: $R2"; exit 1; }

# job-private workspace keyed by the SLURM job id -> two concurrent jobs on the
# same sample can never share a path. CoverM's intermediate BAMs go here (TMPDIR)
# and are discarded on exit.
TMP="${WORK}/tmp_coverm/${s}.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
rm -rf "$TMP"; mkdir -p "${TMP}/out"
export TMPDIR="$TMP"
trap 'rm -rf "$TMP"' EXIT

THREADS="${SLURM_CPUS_PER_TASK:-32}"

echo "CoverM contig (minimap2-sr) map+quant ... $(date +'%H:%M:%S')"
coverm contig \
    -1 "$R1" \
    -2 "$R2" \
    --reference "$CATALOG" \
    --mapper minimap2-sr \
    --methods count covered_fraction \
    --min-covered-fraction 0 \
    --threads "$THREADS" \
    --bam-file-cache-directory "$TMP" \
    --output-file "${TMP}/out/${s}.tsv"

ngenes=$(($(wc -l < "${TMP}/out/${s}.tsv") - 1))
echo "=== $s (${comp:-unknown}): ${ngenes} catalog genes reported ==="

# atomic publish: only a fully-successful run reaches ${OUT}/${s}.tsv
mv "${TMP}/out/${s}.tsv" "${OUT}/${s}.tsv.partial"
mv "${OUT}/${s}.tsv.partial" "${OUT}/${s}.tsv"

touch "${DONE}/${s}.done"
echo "Done $s at $(date +'%F %H:%M:%S')"
