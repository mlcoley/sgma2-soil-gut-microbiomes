#!/bin/bash
#SBATCH --job-name=sgma2_wol2
#SBATCH --output=logs/wol2_%A_%a.out
#SBATCH --error=logs/wol2_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --array=1-256%20
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Per-sample taxonomic profiling: host-filtered R1+R2 (concatenated, unpaired) are
# aligned with bowtie2 to the Web of Life 2 (WoLr2) index, classified with Woltka at
# phylum-species plus per-genome coverage, and the SAM is deleted.
# One array task per row of $SGMA2_BASE/fastq/sample_list.tsv; set --array=1-N.
# Versions: bowtie2 2.5.5, woltka 0.1.7, WoLr2 index 2022-07-25.
# Output: taxonomy/woltka/<sample>*

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
echo "Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_JOB_ID} starting $(date +'%F %H:%M:%S') on $(hostname)"

set +u
source ~/.bashrc
conda activate metagenomic_env
set -u

THREADS=${SLURM_CPUS_PER_TASK:-16}

# ---- paths ----
BASE="${SGMA2_BASE}"
SAMPLELIST="${BASE}/fastq/sample_list.tsv"
READDIR="${BASE}/host_filter/filtered_fastq"
WORK="${BASE}/taxonomy"
SAMDIR="${WORK}/bowtie2"
WOLDIR="${WORK}/woltka"
DONE_DIR="${WORK}/done"
mkdir -p "$SAMDIR" "$WOLDIR" "$DONE_DIR" logs

# ---- WoL2 database / helpers ----
WOLDB="${SGMA2_DB_DIR}/wol2/bowtie2/WoLr2"
LINEAGE="${SGMA2_DB_DIR}/wol2/taxonomy/lineages.txt"
LENGTHMAP="${SGMA2_DB_DIR}/wol2/genomes/length.map"
COVSCRIPT="${SGMA2_COVERAGE_SCRIPT}"

# ---- pull this task's sample ----
sample_line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLELIST")
if [ -z "$sample_line" ]; then echo "No row ${SLURM_ARRAY_TASK_ID} — past end of list, exiting."; exit 0; fi
samplename=$(echo "$sample_line" | cut -f1)
echo "Sample name: $samplename"

# ---- resume: skip if already done ----
if [ -f "${DONE_DIR}/${samplename}.done" ]; then echo "$samplename already done — exiting."; exit 0; fi

R1="${READDIR}/${samplename}_1.filt.fq.gz"
R2="${READDIR}/${samplename}_2.filt.fq.gz"
[ -f "$R1" ] || { echo "ERROR: R1 missing: $R1"; exit 1; }
[ -f "$R2" ] || { echo "ERROR: R2 missing: $R2"; exit 1; }

SAM="${SAMDIR}/${samplename}.sam"

# ---- 1) Bowtie2 -> WoLr2 (WoL2-recommended parameters) ----
echo "bowtie2 -> WoLr2 ... $(date +'%H:%M:%S')"
zcat "$R1" "$R2" | bowtie2 \
    -p "$THREADS" \
    -x "$WOLDB" \
    -S "$SAM" \
    --seed 42 \
    --very-sensitive -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" \
    --score-min "L,0,-0.05" \
    --no-head --no-unal \
    -q -U -

# ---- 2) Woltka classify (multi-rank tables + coverage) ----
echo "woltka classify ... $(date +'%H:%M:%S')"
woltka classify \
    -i "$SAM" \
    -o "${WOLDIR}/${samplename}" \
    --no-demux \
    --lineage "$LINEAGE" \
    --rank phylum,class,order,family,genus,species,none,free \
    --outcov "${WOLDIR}/${samplename}_coverages"

# ---- 3) Per-genome percent coverage ----
echo "coverage ... $(date +'%H:%M:%S')"
python "$COVSCRIPT" "$LENGTHMAP" \
    "${WOLDIR}/${samplename}_coverages/${samplename}.cov" \
    "${WOLDIR}/${samplename}_processed.cov"

# ---- 4) drop the (huge) SAM ----
rm -f "$SAM"

touch "${DONE_DIR}/${samplename}.done"
echo "Done $samplename at $(date +'%F %H:%M:%S')"
