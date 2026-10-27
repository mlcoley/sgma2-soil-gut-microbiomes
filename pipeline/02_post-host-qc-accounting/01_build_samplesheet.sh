#!/usr/bin/env bash
# Joins host-filtered FASTQs with compartment and extraction metadata by P-number.
# Run once from the repository root, then check the printed summary.
# Output: $SGMA2_BASE/posthost_qc/step1_samplesheet.tsv (sample_id, pnum, field_name, compartment, r1, r2,
#   fastp_json, gdna_conc, pcr_conc)

source "${SGMA2_CONFIG:-config/paths.sh}"

set -euo pipefail

# ---- paths ----
BASE="${SGMA2_BASE}/host_filter"
FILT_DIR="${BASE}/filtered_fastq"          # contains *_1.filt.fq.gz / *_2.filt.fq.gz
FASTP_LOG_DIR="${BASE}/logs"               # contains *.fastp.json
META="${SGMA2_BASE}/646_metag-prep-info.txt"   # the metadata sheet
OUT="${SGMA2_BASE}/posthost_qc/step1_samplesheet.tsv"
mkdir -p "$(dirname "$OUT")"

[ -d "$FILT_DIR" ] || { echo "ERROR: FILT_DIR not found: $FILT_DIR"; exit 1; }
[ -f "$META" ]     || { echo "ERROR: metadata not found: $META"; exit 1; }

# Pre-parse metadata into a lookup keyed by sample P-number (strip CR; handle controls).
#   key=Pxxx -> field_name \t compartment \t gdna \t pcr
META_LOOKUP=$(mktemp)
tr -d '\r' < "$META" | awk -F'\t' '
  NR==1 { next }
  $1=="" { next }
  {
    pnum=""; if (match($1, /P[0-9]{3}/)) pnum=substr($1, RSTART, RLENGTH);
    if (pnum=="") next;
    name=$2; gdna=$5; pcr=$6;
    comp="unknown";
    if      (name ~ /^H[0-9]+-S[0-9]+/) comp="soil";
    else if (name ~ /^H[0-9]+-G[0-9]+/) comp="gut";
    else if (name ~ /^H[0-9]+-L[0-9]+/) comp="dung";
    else if (name ~ /PCRBlank/)         comp="control_blank";
    else if (name ~ /Zymo/)             comp="control_mock";
    else if (name ~ /Blank|blank|Neg|Ext/) comp="control_blank";
    print pnum"\t"name"\t"comp"\t"gdna"\t"pcr;
  }' > "$META_LOOKUP"

echo "metadata rows indexed: $(wc -l < "$META_LOOKUP")"

# Walk the filtered R1 files, derive sample_id + P-number, join metadata.
: > "$OUT"
n=0; nomatch=0
shopt -s nullglob
for r1 in "$FILT_DIR"/*_1.filt.fq.gz; do
    sid=$(basename "$r1" _1.filt.fq.gz)
    r2="${FILT_DIR}/${sid}_2.filt.fq.gz"
    [ -f "$r2" ] || { echo "WARN: no R2 for $sid (skipping)"; continue; }

    pnum=$(echo "$sid" | grep -oE 'P[0-9]{3}' | head -1 || true)
    fastp_json="${FASTP_LOG_DIR}/${sid}.fastp.json"
    [ -f "$fastp_json" ] || fastp_json=""   # leave blank if not found; builder will flag

    if [ -n "$pnum" ]; then
        meta=$(awk -F'\t' -v p="$pnum" '$1==p{print $2"\t"$3"\t"$4"\t"$5; exit}' "$META_LOOKUP")
    else
        meta=""
    fi
    if [ -z "$meta" ]; then
        meta=$'NA\tunknown\tNA\tNA'; nomatch=$((nomatch+1))
        echo "WARN: no metadata match for $sid (pnum='${pnum:-none}')"
    fi
    field=$(echo "$meta" | cut -f1)
    comp=$(echo  "$meta" | cut -f2)
    gdna=$(echo  "$meta" | cut -f3)
    pcr=$(echo   "$meta" | cut -f4)

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sid" "${pnum:-NA}" "$field" "$comp" "$r1" "$r2" "$fastp_json" "$gdna" "$pcr" >> "$OUT"
    n=$((n+1))
done
rm -f "$META_LOOKUP"

echo ""
echo ">> Wrote $n samples to $OUT  (unmatched: $nomatch)"
echo ">> Compartment tally:"
cut -f4 "$OUT" | sort | uniq -c
echo ""
echo ">> First 3 rows:"
head -n 3 "$OUT"
echo ""
echo ">> Check before submitting the array:"
echo "   - tally should read ~163 soil / 43 gut / 43 dung / 3 control_blank / 4 control_mock"
echo "   - unmatched should be 0; fastp_json column should be populated (not empty)"
