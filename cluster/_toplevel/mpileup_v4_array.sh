#!/bin/bash
# =============================================================================
# cluster/_toplevel/mpileup_v4_array.sh
# =============================================================================
# One pileup per pool (array task = line of bamlist_v4.txt), samtools mpileup
# -B -q 10 -Q 20 -d 500 -a against the pipe-form reference. -a writes every
# position, so all 27 files have identical rows and can be pasted side by
# side; -d 500 caps per-pool depth (coverage is ~20-30x). Contig names are
# sanitized (| -> _) inline. First job of submit_v4.sh.
#
# Step: 3, pileup (cluster/README.md)
# Manuscript: Supplementary Methods "Read processing, alignment, variant
#             calling and software versions"
#
# Inputs:  bamlist_v4.txt, locations_v4.txt (27 pools, same order)
#          <pool>/samfile/pip_1021/<pool>.filtered.bam (01_alignment/filter_index_flagstat.sh)
#          reference: REF_FASTA
# Outputs: mpileup_parts_v4/<pool>.mpileup.gz
# =============================================================================
#SBATCH --job-name=mpileup_v4_array
#SBATCH --array=1-27
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH -o logs/mpileup_v4_array.%A_%a.out
#SBATCH -e logs/mpileup_v4_array.%A_%a.err

set -Eeuo pipefail
trap 'echo "[ERROR] line $LINENO status $?"' ERR

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

BASE="${PROJECT_ROOT}"
REF="${REF_FASTA}"
BAMLIST="${BASE}/bamlist_v4.txt"
LOCATIONS="${BASE}/locations_v4.txt"
PARTS_DIR="${BASE}/mpileup_parts_v4"

source "${CONDA_SH}"
conda activate "${ENV_POOLSNP}"
module load samtools

mkdir -p "$PARTS_DIR" logs

# --- Identify this task's pool ----------------------------------------------
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAMLIST")
POOL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LOCATIONS" | tr -d ' \t\r\n')
OUT="${PARTS_DIR}/${POOL}.mpileup.gz"

echo "[INFO] task ${SLURM_ARRAY_TASK_ID}: pool=${POOL}"
echo "[INFO] BAM: ${BAM}"

[[ -s "$BAM" ]] || { echo "[ERROR] missing/empty BAM: $BAM"; exit 1; }
[[ -s "$REF" ]] || { echo "[ERROR] missing REF: $REF";  exit 1; }

# --- Run mpileup ------------------------------------------------------------
samtools mpileup \
  -B \
  -q 10 \
  -Q 20 \
  -d 500 \
  -a \
  -f "$REF" \
  "$BAM" \
  | sed 's/|/_/g' \
  | gzip > "$OUT"

# --- Sanity (uses `read` from process substitution to avoid SIGPIPE under
# --- pipefail; `zcat | head -1 | awk` triggers exit 141 because head closes
# --- its stdin and SIGPIPEs zcat, which pipefail then propagates) -----------
read -r first_line < <(zcat "$OUT")
ncols=$(echo "$first_line" | awk '{print NF}')
chrom=$(echo "$first_line" | cut -f1)
[[ "$ncols" -eq 6 ]] || { echo "[ERROR] expected 6 columns, got $ncols"; exit 1; }
[[ "$chrom" != *"|"* ]] || { echo "[ERROR] pipe-form chrom still present: $chrom"; exit 1; }

echo "[✓] task ${SLURM_ARRAY_TASK_ID} (${POOL}) done: $OUT  (cols=$ncols  chrom=$chrom)"
