#!/bin/bash
#SBATCH --job-name=assembly_setup
#SBATCH --output=logs/assembly_setup_%j.out
#SBATCH --error=logs/assembly_setup_%j.err
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
# Builds assembly conda env (MEGAHIT, metaSPAdes, QUAST, seqkit). Run this once.

set -euo pipefail
set +u
source ~/.bashrc
set -u

ENV=assembly

echo "=== conda version ==="
conda --version

echo "=== create env $ENV $(date +'%F %H:%M:%S') ==="
conda create -y -n "$ENV" -c conda-forge -c bioconda \
    megahit spades quast seqkit

set +u
conda activate "$ENV"
set -u

echo "=== resolved versions ==="
echo -n "megahit:  "; megahit --version
echo -n "spades:   "; spades.py --version
echo -n "quast:    "; quast.py --version 2>&1 | head -1
echo -n "seqkit:   "; seqkit version

echo "assembly setup DONE $(date +'%F %H:%M:%S')"
