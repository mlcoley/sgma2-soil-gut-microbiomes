#!/bin/bash
#SBATCH --job-name=amr_setup
#SBATCH --output=logs/amr_setup_%j.out
#SBATCH --error=logs/amr_setup_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Builds the amr conda env (AMRFinderPlus, seqkit), downloads the AMRFinderPlus
# database and records its version. Run once.

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
set +u; source ~/.bashrc; set -u

BASE="${SGMA2_BASE}"
WORK="${BASE}/assembly"
AMR="${WORK}/genecat/amr"
DB="${AMR}/amrfinder_db"
mkdir -p logs "$AMR" "$DB"

# ---- build env (skip if it already exists) ----
if ! conda env list | grep -qE '(^|/)amr[[:space:]]'; then
  echo "Creating conda env 'amr' (ncbi-amrfinderplus + seqkit) ..."
  conda create -y -n amr -c conda-forge -c bioconda ncbi-amrfinderplus seqkit
else
  echo "conda env 'amr' already exists — skipping create."
fi
conda activate amr

echo "amrfinder version: $(amrfinder --version 2>&1 || true)"

# ---- fetch + pin the reference DB into scratch ----
echo "Downloading AMRFinderPlus DB -> $DB   $(date +'%F %H:%M:%S')"
if command -v amrfinder_update >/dev/null 2>&1; then
  amrfinder_update --force_update -d "$DB"
else
  amrfinder -u -d "$DB" || amrfinder -u   # fallback: update default location
fi

# resolve the 'latest' pointer amrfinder_update lays down and stamp the version
LATEST="${DB}/latest"
[ -e "$LATEST" ] || LATEST=$(ls -d "${DB}"/*/ 2>/dev/null | sort | tail -1)
{
  echo "AMRFinderPlus DB setup  $(date +'%F %H:%M:%S')"
  echo "amrfinder version: $(amrfinder --version 2>&1 || true)"
  echo "DB dir:            $DB"
  echo "DB latest:         $LATEST"
  [ -f "${LATEST%/}/version.txt" ] && echo "DB version:        $(cat "${LATEST%/}/version.txt")"
} | tee "${AMR}/amrfinder_db_version.txt"

echo "Done $(date +'%F %H:%M:%S')"
