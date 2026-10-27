# Site configuration. Edit, then submit every job from the repository root:
#   cd <repo> && mkdir -p logs && sbatch pipeline/<stage>/<script>
# Scripts source this file; override its location with SGMA2_CONFIG.
# SLURM account and email are set in each script's #SBATCH header
# (YOUR_ACCOUNT, YOUR_EMAIL).

# Project scratch root: holds fastq/, host_filter/, posthost_qc/, taxonomy/,
# biobakery/, assembly/ and metadata sheets.
export SGMA2_BASE="${SGMA2_BASE:-/path/to/scratch/sgma2}"

# Shared reference databases (minimap2 indexes, WoL2, MetaPhlAn, eggNOG).
export SGMA2_DB_DIR="${SGMA2_DB_DIR:-/path/to/databases}"

# calculate-coverage.py (per-genome coverage from Woltka --outcov output), copied from qiita-spots/qp-woltka.
export SGMA2_COVERAGE_SCRIPT="${SGMA2_COVERAGE_SCRIPT:-pipeline/03_taxonomy-woltka/calculate-coverage.py}"
