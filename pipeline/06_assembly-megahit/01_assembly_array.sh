#!/bin/bash
#SBATCH --job-name=asm
#SBATCH --output=logs/asm_%A_%a.out
#SBATCH --error=logs/asm_%A_%a.err
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
# Per-sample de novo assembly with MEGAHIT (all compartments), contigs >= 1000 bp,
# then QUAST. Samples are never co-assembled across compartments.
# One array task per row of $SGMA2_BASE/assembly/assembly_manifest.tsv
# (sample, compartment), built from the sample metadata sheet:
#   awk -F'\t' 'NR>1 && $5=="true_sample" && $9>=20000000{print $1"\t"$2}' \
#     posthost_qc/decontam_metadata.tsv | sort -t$'\t' -k2,2 -k1,1 > assembly/assembly_manifest.tsv
# Submit with --array=1-N%20 for N manifest rows. Resources: 32 CPUs, 120 GB.
# Versions: MEGAHIT 1.2.9, QUAST 5.3.0, seqkit 2.13.0.
# Output: assembly/assemblies/<sample>/<sample>.contigs.fa and quast/report.tsv

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
echo "Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} start $(date +'%F %H:%M:%S') on $(hostname)"

set +u
source ~/.bashrc
conda activate assembly
set -u

THREADS=${SLURM_CPUS_PER_TASK:-32}
MEM_GB=$(( ${SLURM_MEM_PER_NODE:-122880} / 1024 ))          # SLURM MB -> GB
MEMBYTES=$(( MEM_GB * 1024 * 1024 * 1024 ))                 # for megahit -m
echo "threads=$THREADS  mem=${MEM_GB}G"

BASE="${SGMA2_BASE}"
WORK="${BASE}/assembly"
READDIR="${BASE}/host_filter/filtered_fastq"
OUT="${WORK}/assemblies"
DONE="${WORK}/done_asm"
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

R1="${READDIR}/${s}_1.filt.fq.gz"
R2="${READDIR}/${s}_2.filt.fq.gz"
[ -f "$R1" ] || { echo "ERROR: R1 missing: $R1"; exit 1; }
[ -f "$R2" ] || { echo "ERROR: R2 missing: $R2"; exit 1; }

# job-private workspace keyed by the SLURM job id -> two concurrent jobs on the
# same sample can never share a path.
TMP="${WORK}/tmp_asm/${s}.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
rm -rf "$TMP"; mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

echo "MEGAHIT ... $(date +'%H:%M:%S')"
# megahit refuses a pre-existing -o dir; TMP/mh does not exist yet.
megahit -1 "$R1" -2 "$R2" -o "${TMP}/mh" -t "$THREADS" -m "$MEMBYTES" \
        --min-contig-len 1000 --out-prefix "$s"

# length-filter >=1000 bp (megahit already did via --min-contig-len; harmless
# re-filter pins the final contigs to the same explicit rule)
mkdir -p "${TMP}/out"
CONTIGS="${TMP}/out/${s}.contigs.fa"
seqkit seq --min-len 1000 "${TMP}/mh/${s}.contigs.fa" > "$CONTIGS"

echo "QUAST ... $(date +'%H:%M:%S')"
quast.py "$CONTIGS" -o "${TMP}/out/quast" -t "$THREADS" --min-contig 1000

echo "=== QUAST summary ($s, $comp) ==="
cat "${TMP}/out/quast/report.tsv" || true

# atomic publish: only a fully-successful run reaches ${OUT}/${s}
rm -rf "${OUT}/${s}.partial" "${OUT}/${s}"
mv "${TMP}/out" "${OUT}/${s}.partial"
mv "${OUT}/${s}.partial" "${OUT}/${s}"

touch "${DONE}/${s}.done"
echo "Done $s at $(date +'%F %H:%M:%S')"
