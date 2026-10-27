#!/bin/bash
#SBATCH --job-name=sgma2_s1_agg
#SBATCH --output=logs/s1_agg_%j.out
#SBATCH --error=logs/s1_agg_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A YOUR_ACCOUNT
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL
# Run once after the QC array finishes: MultiQC over fastp and FastQC outputs, then
# the read-accounting table (raw, post-fastp, post-host counts).

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail
set +u
source ~/.bashrc
set -u
conda activate metagenomic_env

WORK="${SGMA2_BASE}/posthost_qc"
SHEET="${WORK}/step1_samplesheet.tsv"
SEQKIT_DIR="${WORK}/qc/seqkit"
FASTQC_DIR="${WORK}/qc/fastqc"
FASTP_LOG_DIR="${SGMA2_BASE}/host_filter/logs"   # *.fastp.json live here
MQC_DIR="${WORK}/qc/multiqc"
ACCT="${WORK}/sgma2_step1_read_accounting.tsv"
BUILDER="pipeline/02_post-host-qc-accounting/build_accounting_table.py"
mkdir -p "$MQC_DIR" logs

# completeness check
n_sheet=$(wc -l < "$SHEET")
n_seqkit=$(ls "${SEQKIT_DIR}"/*.seqkit.tsv 2>/dev/null | wc -l)
echo "samplesheet rows: $n_sheet ; seqkit outputs: $n_seqkit"
[ "$n_seqkit" -lt "$n_sheet" ] && echo "WARNING: fewer seqkit outputs than samples — some array tasks may have failed. Check logs/ before trusting the table."

# ---- 1) MultiQC (auto-detects fastp JSON + FastQC) ----
echo "MultiQC..."
multiqc -f -o "$MQC_DIR" -n sgma2_step1_multiqc "$FASTQC_DIR" "$FASTP_LOG_DIR"

# ---- 2) Master read-accounting table ----
echo "Building read-accounting table..."
python3 "$BUILDER" --sheet "$SHEET" --seqkit "$SEQKIT_DIR" --logdir "$FASTP_LOG_DIR" --out "$ACCT"

echo ""
echo ">> Outputs:"
echo "   MultiQC report:   ${MQC_DIR}/sgma2_step1_multiqc.html"
echo "   Accounting table: ${ACCT}"
echo "Done $(date)"
