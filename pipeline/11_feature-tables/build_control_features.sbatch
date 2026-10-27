#!/bin/bash
#SBATCH --job-name=ctrl_features
#SBATCH --output=logs/ctrl_features_%j.out
#SBATCH --error=logs/ctrl_features_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Runs build_control_features.py after the control libraries have been mapped to the
# gene catalog (08_abundance-coverm/01_coverm_array.sh with MANIFEST=control_manifest.tsv).
# Needs about 48 GB: the gene-to-KO map is held in memory.
# Output: assembly/genecat/controls/

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
echo "start $(date +'%F %H:%M:%S') on $(hostname)"

set +u
source ~/.bashrc
conda activate genecat
set -u

mkdir -p logs
python pipeline/11_feature-tables/build_control_features.py

echo "done $(date +'%F %H:%M:%S')"
