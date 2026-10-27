#!/bin/bash
#SBATCH --job-name=microbecensus
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=1-194%20
#SBATCH --output=logs/microbecensus_%A_%a.out
#SBATCH --error=logs/microbecensus_%A_%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# MicrobeCensus average genome size and genome equivalents per sample, from the
# host-filtered reads (2,000,000 reads subsampled per sample). The merge step builds
# feature_tables/microbecensus_ge.tsv for the analysis.
# One array task per line of $SGMA2_BASE/microbecensus/analysis_samples.txt (sample
# ids of the 194-sample analysis cohort); submit with --array=1-194%20.
# Environment: conda create -n microbecensus2 -c bioconda -c conda-forge microbecensus
# Output: microbecensus/per_sample/<sample>.mc.txt

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
set +u; source ~/.bashrc; set -u
conda activate microbecensus2

BASE=${SGMA2_BASE}
READS=$BASE/host_filter/filtered_fastq
OUT=$BASE/microbecensus
MANIFEST=${MANIFEST:-$OUT/analysis_samples.txt}
mkdir -p "$OUT/per_sample" "$OUT/done" logs

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
[ -z "$SAMPLE" ] && { echo "no sample at row $SLURM_ARRAY_TASK_ID"; exit 1; }
echo "Sample: $SAMPLE"
[ -f "$OUT/done/${SAMPLE}.done" ] && { echo "already done"; exit 0; }

R1=$READS/${SAMPLE}_1.filt.fq.gz
R2=$READS/${SAMPLE}_2.filt.fq.gz
TMP=$(mktemp -d "$OUT/tmp_${SAMPLE}_XXXX")
trap 'rm -rf "$TMP"' EXIT

# MicrobeCensus reads both mates; -n subsamples reads (default 2M is plenty for AGS)
run_microbe_census.py -t "$SLURM_CPUS_PER_TASK" -n 2000000 \
  "$R1","$R2" "$TMP/${SAMPLE}.mc.txt"

mv "$TMP/${SAMPLE}.mc.txt" "$OUT/per_sample/${SAMPLE}.mc.txt"
touch "$OUT/done/${SAMPLE}.done"
echo "done $SAMPLE"
