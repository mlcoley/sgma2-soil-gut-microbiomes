#!/bin/bash
#SBATCH --job-name=mpa
#SBATCH --output=logs/mpa_%A_%a.out
#SBATCH --error=logs/mpa_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --array=1-256%20
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Per-sample MetaPhlAn 4.1.2 profiling of host-filtered reads (R1 and R2 passed
# together; pairing is ignored). One array task per row of
# $SGMA2_BASE/fastq/sample_list.tsv; set --array=1-N.
# Database: mpa_vJun23_CHOCOPhlAnSGB_202307 under $SGMA2_DB_DIR/metaphlan.
# Versions: MetaPhlAn 4.1.2, bowtie2 2.5.5.
# Output: biobakery/metaphlan/<sample>_metaphlan.txt and <sample>.bowtie2.bz2

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
echo "Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} start $(date +'%F %H:%M:%S') on $(hostname)"

set +u
source ~/.bashrc
conda activate biobakery
set -u

THREADS=${SLURM_CPUS_PER_TASK:-16}

BASE="${SGMA2_BASE}"
SAMPLELIST="${BASE}/fastq/sample_list.tsv"
READDIR="${BASE}/host_filter/filtered_fastq"
WORK="${BASE}/biobakery"
OUT="${WORK}/metaphlan"
DONE="${WORK}/done_mpa"
mkdir -p "$OUT" "$DONE" logs

MPA_DB_DIR="${SGMA2_DB_DIR}/metaphlan"
MPA_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202307"

# ---- this task's sample (sample_list line N -> P###) ----
s=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLELIST" | cut -f1)
[ -z "$s" ] && { echo "No row ${SLURM_ARRAY_TASK_ID} — past end, exit."; exit 0; }
echo "Sample name: $s"
[ -f "${DONE}/${s}.done" ] && { echo "$s already done — exit."; exit 0; }

R1="${READDIR}/${s}_1.filt.fq.gz"
R2="${READDIR}/${s}_2.filt.fq.gz"
[ -f "$R1" ] || { echo "ERROR: R1 missing: $R1"; exit 1; }
[ -f "$R2" ] || { echo "ERROR: R2 missing: $R2"; exit 1; }

# clean any partial intermediates from a prior failed run (metaphlan refuses an
# existing bowtie2out); done-marker gates legitimate resumes.
rm -f "${OUT}/${s}.bowtie2.bz2" "${OUT}/${s}_metaphlan.txt"

echo "metaphlan ... $(date +'%H:%M:%S')"
metaphlan "${R1},${R2}" \
    --input_type fastq \
    --bowtie2db "$MPA_DB_DIR" --index "$MPA_INDEX" \
    --bowtie2out "${OUT}/${s}.bowtie2.bz2" \
    --nproc "$THREADS" \
    --unclassified_estimation \
    --sample_id "$s" \
    -o "${OUT}/${s}_metaphlan.txt"

touch "${DONE}/${s}.done"
echo "Done $s at $(date +'%F %H:%M:%S')"
