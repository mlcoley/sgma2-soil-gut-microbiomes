#!/bin/bash
#SBATCH --job-name=biobakery_setup
#SBATCH --output=logs/biobakery_setup_%j.out
#SBATCH --error=logs/biobakery_setup_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Builds a biobakery conda env (MetaPhlAn 4.1.2) and checks the MetaPhlAn index
# under $SGMA2_DB_DIR/metaphlan. Run once

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
set +u
source ~/.bashrc
set -u

LABDB="${SGMA2_DB_DIR}"
ENV=biobakery

echo "=== conda version ==="
conda --version

echo "=== create env $ENV (metaphlan=4.1.2) $(date +'%F %H:%M:%S') ==="
conda create -y -n "$ENV" -c conda-forge -c bioconda metaphlan=4.1.2

set +u
conda activate "$ENV"
set -u

echo "=== resolved versions ==="
echo -n "metaphlan: "; metaphlan --version
echo -n "bowtie2:   "; bowtie2 --version | head -1
echo -n "diamond:   "; diamond --version 2>/dev/null || echo "(n/a)"

echo "=== MetaPhlAn index visible (no download) ==="
ls "$LABDB/metaphlan/" | grep mpa_vJun23_CHOCOPhlAnSGB_202307

echo "biobakery setup DONE $(date +'%F %H:%M:%S')"
