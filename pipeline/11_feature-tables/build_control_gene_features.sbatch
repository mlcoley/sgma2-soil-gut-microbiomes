#!/bin/bash
#SBATCH --job-name=ctrl_genes
#SBATCH --output=logs/ctrl_genes_%j.out
#SBATCH --error=logs/ctrl_genes_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Runs build_control_gene_features.py on the control CoverM outputs.
# Output: assembly/genecat/controls/

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail

mkdir -p logs

set +u
source ~/.bashrc
conda activate genecat
set -u

echo "host: $(hostname)   start: $(date)"
python pipeline/11_feature-tables/build_control_gene_features.py
echo "end: $(date)"
