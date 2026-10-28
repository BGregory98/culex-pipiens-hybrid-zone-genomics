#!/bin/bash
# =============================================================================
# cluster/04_baypass/vcf_to_allele_counts.sh
# =============================================================================
# Converts the PoolSNP VCF to a per-SNP table of per-pool reference and
# alternate read counts, the input of prepare_baypass_input_v4.R. Run
# interactively from the project root with the bcftools module loaded.
#
# Step: 6, BayPass (cluster/README.md)
# Inputs:  PoolSNP_out_v4/poolsnp.biallelic.nomiss.vcf.gz (step 5)
# Outputs: allele_counts_v4.tsv (also read by the R stage; deposit
#          data/allele_counts_v4.tsv)
#
# Output format (TSV, no header):
#   col 1     CHROM (chromosome name canonicalized to underscore form
#             ENA_OZ00xxxx_OZ00xxxx.1, matching the downstream pipeline)
#   col 2     POS
#   col 3     REF allele
#   col 4     ALT allele
#   cols 5-N  per-pool "RD:AD" strings (one column per sample in VCF order)
#
# PoolSNP-specific note: in PoolSNP VCFs, AD is the alt count alone (Number=1)
# and RD is the ref count -- this is non-standard but well-documented. We pull
# both fields directly via bcftools query.
#
# Usage:
#   ./vcf_to_allele_counts.sh <input.vcf.gz> <output.tsv>
#
# Run as:
#   bash scripts/04_baypass/vcf_to_allele_counts.sh \
#     PoolSNP_out_v4/poolsnp.biallelic.nomiss.vcf.gz \
#     allele_counts_v4.tsv
#
# Requires: bcftools (load via `module load bcftools` or equivalent).
# =============================================================================

set -Eeuo pipefail

VCF=${1:?Usage: vcf_to_allele_counts.sh <input.vcf.gz> <output.tsv>}
OUT=${2:?Usage: vcf_to_allele_counts.sh <input.vcf.gz> <output.tsv>}

# ---- Sanity checks --------------------------------------------------------

if [[ ! -s "$VCF" ]]; then
  echo "ERROR: VCF not found or empty: $VCF" >&2
  exit 1
fi

if ! command -v bcftools >/dev/null 2>&1; then
  echo "ERROR: bcftools not in PATH. Try: module load bcftools" >&2
  exit 1
fi

# ---- Confirm AD/RD are in the FORMAT header (PoolSNP-style) --------------

if ! bcftools view -h "$VCF" | grep -q 'ID=AD,'; then
  echo "ERROR: VCF FORMAT header does not contain AD field. Is this a PoolSNP VCF?" >&2
  exit 1
fi
if ! bcftools view -h "$VCF" | grep -q 'ID=RD,'; then
  echo "ERROR: VCF FORMAT header does not contain RD field. Is this a PoolSNP VCF?" >&2
  exit 1
fi

# ---- Print sample order so the user can confirm it matches samples.txt ---

echo "[INFO] Sample order in VCF (preserve this in samples.txt):"
bcftools query -l "$VCF" | nl

# ---- Extract counts; canonicalize chromosome names to underscore form ----
#
# The awk step replaces every "|" in column 1 with "_" so that the chromosome
# names match what the downstream R scripts expect. Pipes appear only in
# chromosome names (BAM-style "ENA|OZ00xxxx|OZ00xxxx.1"); cols 2-N never
# contain "|".

echo "[INFO] Extracting counts to: $OUT"

bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT[\t%RD:%AD]\n' \
  "$VCF" \
  | awk 'BEGIN{OFS="\t"} {gsub(/\|/, "_", $1); print}' \
  > "$OUT"

# ---- Summary ---------------------------------------------------------------

NLINES=$(wc -l < "$OUT")
NCOLS=$(awk -F'\t' '{print NF; exit}' "$OUT")
NPOOLS=$((NCOLS - 4))

echo "[✓] Wrote $OUT"
echo "    Variants: $NLINES"
echo "    Columns:  $NCOLS  ($NPOOLS pools + 4 metadata)"
echo "    Chromosomes seen:"
awk -F'\t' '{print "      ", $1}' "$OUT" | sort -u
