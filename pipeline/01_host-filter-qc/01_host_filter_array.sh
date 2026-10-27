#!/bin/bash
#SBATCH --job-name=sgma2_hostfilter
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --array=1-256%20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --mem=48G
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Per-sample host depletion: fastp adapter/length trim, then minimap2 against the
# cattle (Bos taurus + indicus), GRCh38 and T2T-Y indexes, keeping read pairs that
# are unmapped at each stage. pairfilt.awk removes zero-length reads so pairs stay
# synchronised. Reads row N of $SGMA2_BASE/fastq/sample_list.tsv
# (sample, R1 path, R2 path); set --array=1-N for N samples.
# Output: host_filter/filtered_fastq/<sample>_{1,2}.filt.fq.gz
# SRA_HUMAN_ONLY=true (default false) gives the SRA deposit: no adapter, quality or
# length trimming and no cattle step; reads are interleaved by fastp, then filtered
# against GRCh38 and T2T-Y only. Output: host_filter_human_only/filtered_fastq/
# (same file names). Submit with: SRA_HUMAN_ONLY=true sbatch ...

source "${SGMA2_CONFIG:-config/paths.sh}"

echo "Starting $SLURM_JOBID at $(date +"%H:%M.%S")"

source ~/.bashrc

conda activate metagenomic_env

# Stop at the first failed command or pipeline stage so a failed sample never reaches the done marker.
set -eo pipefail
trap 'echo "ERROR: sample ${samplename:-unknown} failed at line $LINENO" >&2' ERR

# Parse sample list
sample_line=$(head -n ${SLURM_ARRAY_TASK_ID} ${SGMA2_BASE}/fastq/sample_list.tsv | tail -n 1)

echo "Full sample line: $sample_line"

samplename=$(echo "$sample_line" | cut -f1)
echo "Sample name: $samplename"

SRA_HUMAN_ONLY="${SRA_HUMAN_ONLY:-false}"
base_dir="${SGMA2_BASE}/host_filter"
[ "$SRA_HUMAN_ONLY" = true ] && base_dir="${SGMA2_BASE}/host_filter_human_only"
PAIRFILT="pipeline/01_host-filter-qc/pairfilt.awk"

done_dir="${base_dir}/done"
todo_dir="${base_dir}/todo"
prog_dir="${base_dir}/prog"
mkdir -p "$done_dir" "$todo_dir" "$prog_dir"

if [ -f "${done_dir}/${samplename}.txt" ]; then
    echo "Sample $samplename already done, exiting..."
    exit 0
fi

if [ -f "${todo_dir}/${samplename}.txt" ]; then
    mv "${todo_dir}/${samplename}.txt" "${prog_dir}/${samplename}.txt"
elif [ ! -f "${prog_dir}/${samplename}.txt" ]; then
    echo "x" > "${prog_dir}/${samplename}.txt"
fi

forward=$(echo "$sample_line" | cut -f2)
echo "R1 file: $forward"

reverse=$(echo "$sample_line" | cut -f3)
echo "R2 file: $reverse"

# Sanity checks
ls $forward
ls $reverse
file $forward
file $reverse

use_combined=1
decompress=1

LOG_DIR="${base_dir}/logs"
OUT_DIR="${base_dir}/filtered_fastq"

DB_GRCH38="${SGMA2_DB_DIR}/human/grch38-db-may2024.mmi"
DB_T2TY="${SGMA2_DB_DIR}/human/t2ty-db.mmi"
DB_COW="${SGMA2_DB_DIR}/cow/taurus-indicus-db.mmi"

TAURUS_DB="${SGMA2_DB_DIR}/cow/taurus-db.mmi"
INDICUS_DB="${SGMA2_DB_DIR}/cow/indicus-db.mmi"

case "${SGMA2_DB_DIR:-}" in
    ""|/path/to/*) echo "ERROR: SGMA2_DB_DIR is not set; edit config/paths.sh or export it" >&2; exit 1 ;;
esac
needed="$DB_GRCH38 $DB_T2TY"
[ "$SRA_HUMAN_ONLY" = true ] || needed="$needed $DB_COW"
for db in $needed; do
    [ -f "$db" ] || { echo "ERROR: minimap2 index not found: $db" >&2; exit 1; }
done

mkdir -p $LOG_DIR $OUT_DIR

t0=$(date +"%s")

tmp_dir="${base_dir}/tmp"
mkdir -p "$tmp_dir"

fwd_tmp="${tmp_dir}/${samplename}_fwd.fastq"
rev_tmp="${tmp_dir}/${samplename}_rev.fastq"

# --- Decompress ---
if [ "$decompress" -eq 1 ]; then

    [ ! -f "$fwd_tmp" ] && unpigz -c -d -p 16 "$forward" > "$fwd_tmp"
    [ ! -f "$rev_tmp" ] && unpigz -c -d -p 16 "$reverse" > "$rev_tmp"

else
    fwd_tmp="$forward"
    rev_tmp="$reverse"
fi

echo "Decompress done at $(date +"%H:%M.%S")"

if [ "$SRA_HUMAN_ONLY" = false ]; then
    # --- Run fastp ----
    echo "Running fastp..."
    fastp -i "$fwd_tmp" -I "$rev_tmp" \
        --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -l 100 -w 2 -z 1 \
        -j "${LOG_DIR}/${samplename}.fastp.json" \
        -h "${LOG_DIR}/${samplename}.fastp.html" \
        --stdout > "${OUT_DIR}/${samplename}_interleaved.fq"

    echo "Done: fastp at $(date +"%H:%M.%S")"

    if [ "$use_combined" -eq 1 ]; then
        # --- Filter on Cow (taurus + indicus) ----
        echo "Running minimap2 with combined cow db..."
        #  Map to combined cow DB, filter for unmapped paired reads (-f 12), and output as interleaved fastq
        minimap2 -t 16 -ax sr "$DB_COW" "${OUT_DIR}/${samplename}_interleaved.fq" \
            | samtools sort -n - \
            | samtools fastq -@ 16 -f 12 -F 256 -0 /dev/null - \
            | awk -f "${PAIRFILT}" > "${OUT_DIR}/${samplename}_cow_intl.fq"

    else
        echo "Running minimap2 with separate cow dbs..."
        minimap2 -t 16 -ax sr "$TAURUS_DB" "${OUT_DIR}/${samplename}_interleaved.fq" \
            | samtools sort -n - \
            | samtools fastq -@ 16 -f 12 -F 256 -0 /dev/null - \
            | awk -f "${PAIRFILT}" > "${OUT_DIR}/${samplename}_cow1_intl.fq"

        minimap2 -t 16 -ax sr "$INDICUS_DB" "${OUT_DIR}/${samplename}_cow1_intl.fq" \
            | samtools sort -n - \
            | samtools fastq -@ 16 -f 12 -F 256 -0 /dev/null - \
            | awk -f "${PAIRFILT}" > "${OUT_DIR}/${samplename}_cow_intl.fq"

        rm "${OUT_DIR}/${samplename}_cow1_intl.fq"
    fi

    echo "Done: cow filter at $(date +"%H:%M.%S")"

    current_file="${OUT_DIR}/${samplename}_cow_intl.fq"
else
    echo "Running fastp (interleave only, no trimming or filtering)..."
    fastp -i "$fwd_tmp" -I "$rev_tmp" \
        --disable_quality_filtering \
        --disable_length_filtering \
        --disable_adapter_trimming \
        -j "${LOG_DIR}/${samplename}.fastp.json" \
        -h "${LOG_DIR}/${samplename}.fastp.html" \
        --stdout > "${OUT_DIR}/${samplename}_interleaved.fq"

    current_file="${OUT_DIR}/${samplename}_interleaved.fq"
fi

# --- Filter on GRCh38 ----
echo "Running minimap2 with GRCh38..."
#  Map to GRCh38, filter for unmapped paired reads (-f 12), and output as interleaved fastq
minimap2 -t 16 -ax sr "$DB_GRCH38" "$current_file" \
    | samtools sort -n - \
    | samtools fastq -@ 16 -f 12 -F 256 -0 /dev/null - \
    | awk -f "${PAIRFILT}" > "${OUT_DIR}/${samplename}_grch38_intl.fq"

echo "Done: GRCh38 filter at $(date +"%H:%M.%S")"

# --- Filter on T2TY (no PhiX) ----
echo "Running minimap with T2TY..."

# Map remaining reads to T2TY + PhiX, filter for unmapped paired reads (-f 12), and write final compressed files
minimap2 -t 16 -a -x sr "$DB_T2TY" "${OUT_DIR}/${samplename}_grch38_intl.fq" \
    | samtools sort -n - \
    | samtools fastq -@ 16 -f 12 -F 256 -0 /dev/null \
	    -1 "${OUT_DIR}/${samplename}_1.filt.fq.gz" \
 	    -2 "${OUT_DIR}/${samplename}_2.filt.fq.gz" -

echo "Done: T2T-Y filter at $(date +"%H:%M.%S")"

# Remove the intermediate interleaved file
rm -f "${OUT_DIR}/${samplename}_interleaved.fq" "${OUT_DIR}/${samplename}_cow_intl.fq" "${OUT_DIR}/${samplename}_grch38_intl.fq"

if [ "$decompress" -eq 1 ]; then
    rm "$fwd_tmp" "$rev_tmp"
fi

echo "$samplename pipeline complete at $(date +"%H%M")"

mv "${prog_dir}/${samplename}.txt" "${done_dir}/${samplename}.txt"

t1=$(date +"%s")

dur=0
dur=$(( t1 - t0 ))
echo "Full pipeline took $dur seconds"
