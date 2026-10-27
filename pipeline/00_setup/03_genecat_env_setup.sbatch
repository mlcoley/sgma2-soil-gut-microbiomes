#!/bin/bash
#SBATCH --job-name=genecat_setup
#SBATCH --output=logs/genecat_setup_%j.out
#SBATCH --error=logs/genecat_setup_%j.err
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
# Builds the genecat conda env (prodigal, mmseqs2, coverm, eggnog-mapper, seqkit).
# The eggNOG database is not downloaded here. Run once.

set -euo pipefail
set +u
source ~/.bashrc
set -u

ENV=genecat

echo "=== conda version ==="
conda --version

echo "=== create env $ENV $(date +'%F %H:%M:%S') ==="
conda create -y -n "$ENV" -c conda-forge -c bioconda \
    prodigal mmseqs2 coverm eggnog-mapper seqkit

set +u
conda activate "$ENV"
set -u

echo "=== resolved versions ==="
echo -n "prodigal:      "; prodigal -v 2>&1 | head -1
echo -n "mmseqs2:       "; mmseqs version 2>&1 | head -1
echo -n "coverm:        "; coverm --version 2>&1 | head -1
echo -n "eggnog-mapper: "; emapper.py --version 2>&1 | head -1
echo -n "seqkit:        "; seqkit version

echo "genecat setup DONE $(date +'%F %H:%M:%S')"
