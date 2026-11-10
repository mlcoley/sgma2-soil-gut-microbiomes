#!/bin/bash
#SBATCH --job-name=sgma2_s1_qc
#SBATCH --output=logs/s1_qc_%A_%a.out
#SBATCH --error=logs/s1_qc_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --array=1-256%20
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Per-sample read statistics (seqkit) and FastQC on host-filtered reads.
# One array task per row of step1_samplesheet.tsv; set --array=1-N where N is the
# number of rows

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
echo "Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} starting $(date +"%H:%M.%S")"

set +u
source ~/.bashrc
conda activate metagenomic_env
set -u

THREADS=${SLURM_CPUS_PER_TASK:-4}

# ---- paths (edit if needed) ----
WORK="${SGMA2_BASE}/posthost_qc"
SHEET="${WORK}/step1_samplesheet.tsv"
SEQKIT_DIR="${WORK}/qc/seqkit"
FASTQC_DIR="${WORK}/qc/fastqc"
DONE_DIR="${WORK}/done"
mkdir -p "$SEQKIT_DIR" "$FASTQC_DIR" "$DONE_DIR" logs

# ---- pull this task's row ----
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SHEET")
if [ -z "$line" ]; then echo "No row ${SLURM_ARRAY_TASK_ID} (past end of sheet) — exiting."; exit 0; fi

sid=$(echo  "$line" | cut -f1)
comp=$(echo "$line" | cut -f4)
r1=$(echo   "$line" | cut -f5)
r2=$(echo   "$line" | cut -f6)
echo "Sample: $sid   compartment: $comp"

# ---- resume: skip if already done ----
if [ -f "${DONE_DIR}/${sid}.done" ]; then echo "$sid already done — exiting."; exit 0; fi

# sanity
[ -f "$r1" ] || { echo "ERROR: R1 missing: $r1"; exit 1; }
[ -f "$r2" ] || { echo "ERROR: R2 missing: $r2"; exit 1; }

# ---- 1) read/base stats on the host-filtered reads ----
echo "seqkit stats..."
seqkit stats -T -a -j "$THREADS" "$r1" "$r2" > "${SEQKIT_DIR}/${sid}.seqkit.tsv"

# ---- 2) FastQC ----
echo "FastQC..."
fastqc -t 2 -o "$FASTQC_DIR" "$r1" "$r2"

touch "${DONE_DIR}/${sid}.done"
echo "Done $sid at $(date)"
