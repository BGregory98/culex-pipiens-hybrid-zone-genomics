#!/bin/bash
# =============================================================================
# cluster/02_mpileup/sanitize_headers.sh
# =============================================================================
# Part 1 writes a copy of the reference with / | , : replaced by _ in contig
# names (ENA|OZ004311|OZ004311.1 -> ENA_OZ004311_OZ004311.1) and indexes it.
# That FASTA was later moved into reference_genomes/cxpip_genome/ and is the
# reference given to PoolSNP (REF_FASTA_SANITIZED). Part 2 applied the same
# renaming to the combined pileup of an earlier (pre-v4) pileup run, which
# the published analyses do not use; the v4 pileup sanitizes inline
# (_toplevel/mpileup_v4_array.sh). Submitted through run_sanitize.sh.
#
# Step: 4, variant calling (reference preparation; cluster/README.md)
#
# Inputs:  reference_genomes/cxpip_genome.fasta
#          all_pools.mpileup.gz (pre-v4 pileup; part 2 only)
# Outputs: reference_genomes/cxpip_genome.sanitized.fa (+ .fai)
#          all_pools.sanitized.mpileup.gz (+ .tbi; not used downstream)
# =============================================================================
set -euo pipefail

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

module load samtools/1.16 2>/dev/null || true

BASE="${PROJECT_ROOT}"

echo "[INFO] Starting header sanitization at $(date)"
echo "[INFO] Working in: $BASE"

# 1. Sanitize FASTA headers (replace / | , : with _)
awk 'BEGIN{OFS=""} /^>/{gsub(/[\/|,:]/,"_"); print; next} {print}' \
  "$BASE/reference_genomes/cxpip_genome.fasta" > \
  "$BASE/reference_genomes/cxpip_genome.sanitized.fa"

samtools faidx "$BASE/reference_genomes/cxpip_genome.sanitized.fa"

# 2. Sanitize contig names in mpileup (first column only)
zcat "$BASE/all_pools.mpileup.gz" \
| awk 'BEGIN{OFS="\t"}{gsub(/[\/|,:]/,"_",$1); print}' \
| bgzip -@ 4 > "$BASE/all_pools.sanitized.mpileup.gz"

tabix -s 1 -b 2 -e 2 "$BASE/all_pools.sanitized.mpileup.gz"

# 3. Verify results
grep '^>' "$BASE/reference_genomes/cxpip_genome.sanitized.fa" | head
tabix -l "$BASE/all_pools.sanitized.mpileup.gz"

echo "[✓] Sanitization finished successfully at $(date)"

