#!/bin/bash
#SBATCH --job-name=sgma2_merge
#SBATCH --output=logs/merge_%j.out
#SBATCH --error=logs/merge_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Merges per-sample Woltka tables into cohort matrices. Run after the Woltka array
# finishes. Output: taxonomy/merged/

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
set +u
source ~/.bashrc
conda activate metagenomic_env
set -u

python pipeline/03_taxonomy-woltka/merge_woltka_tables.py \
    --woltka-dir ${SGMA2_BASE}/taxonomy/woltka \
    --read-accounting ${SGMA2_BASE}/posthost_qc/sgma2_step1_read_accounting.tsv \
    --out ${SGMA2_BASE}/taxonomy/merged
echo "merge done $(date +'%F %H:%M:%S')"
