#!/bin/bash
# =============================================================================
# cluster/_toplevel/mpileup_v4_merge.sh
# =============================================================================
# Pastes the 27 per-pool pileups into one multi-pool pileup after checking
# that every part exists and has the same number of lines. Columns: chrom,
# pos, ref (from pool 1), then depth/bases/quals for each pool in
# locations_v4.txt order (3 + 3 x 27 = 84). Second job of submit_v4.sh.
#
# Step: 3, pileup (cluster/README.md)
# Manuscript: Supplementary Methods "Read processing, alignment, variant
#             calling and software versions"
#
# Inputs:  mpileup_parts_v4/<pool>.mpileup.gz (mpileup_v4_array.sh)
#          locations_v4.txt
# Outputs: pools27_v4.sanitized.mpileup.gz (read by PoolSNP_v4.sh and the
#          06_depth/glmm callable-depth steps)
# =============================================================================
#SBATCH --job-name=mpileup_v4_merge
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH --time=0-22:00:00
#SBATCH -o logs/mpileup_v4_merge.%j.out
#SBATCH -e logs/mpileup_v4_merge.%j.err

set -Eeuo pipefail
trap 'echo "[ERROR] line $LINENO status $?"' ERR

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

BASE="${PROJECT_ROOT}"
LOCATIONS="${BASE}/locations_v4.txt"
PARTS_DIR="${BASE}/mpileup_parts_v4"
OUT="${BASE}/pools27_v4.sanitized.mpileup.gz"

mkdir -p logs

# --- Confirm all 27 part files exist and are non-empty ----------------------
mapfile -t pools < <(sed 's/[[:space:]]//g' "$LOCATIONS")
echo "[INFO] merging ${#pools[@]} pools"
for pool in "${pools[@]}"; do
  f="${PARTS_DIR}/${pool}.mpileup.gz"
  [[ -s "$f" ]] || { echo "[ERROR] missing part: $f"; exit 1; }
done
echo "[INFO] all part files present"

# --- Confirm all files have the same line count ----------------------------
nlines_ref=$(zcat "${PARTS_DIR}/${pools[0]}.mpileup.gz" | wc -l)
echo "[INFO] ${pools[0]}: $nlines_ref lines"
for pool in "${pools[@]:1}"; do
  n=$(zcat "${PARTS_DIR}/${pool}.mpileup.gz" | wc -l)
  [[ "$n" -eq "$nlines_ref" ]] || {
    echo "[ERROR] line count mismatch: ${pool} has $n lines, expected $nlines_ref"
    exit 1
  }
done
echo "[INFO] all files have $nlines_ref lines — safe to paste"

# --- Build paste command dynamically ----------------------------------------
# Pool 1: all 6 columns (chrom, pos, ref, depth, bases, quals)
# Pools 2-27: columns 4-6 only (depth, bases, quals)
paste_args=("<(zcat \"${PARTS_DIR}/${pools[0]}.mpileup.gz\")")
for pool in "${pools[@]:1}"; do
  paste_args+=("<(zcat \"${PARTS_DIR}/${pool}.mpileup.gz\" | cut -f4-6)")
done

echo "[INFO] pasting ${#pools[@]} files..."
eval "paste ${paste_args[*]}" | gzip > "$OUT"

# --- Sanity (read first line via process substitution to avoid SIGPIPE
# --- under pipefail; see mpileup_v4_array.sh) -------------------------------
read -r first_line < <(zcat "$OUT")
ncols=$(echo "$first_line" | awk '{print NF}')
chrom=$(echo "$first_line" | cut -f1)
echo "[INFO] output columns: $ncols (expected 84)"
echo "[INFO] first chrom: $chrom"
[[ "$ncols" -eq 84 ]] || { echo "[ERROR] column count wrong"; exit 1; }
[[ "$chrom" != *"|"* ]] || { echo "[ERROR] pipe-form chrom in output"; exit 1; }

echo "[✓] merge done: $OUT"
