#!/bin/bash
# =============================================================================
# cluster/_toplevel/PoolSNP_v4.sh
# =============================================================================
# Calls SNPs across the 27 pools (24 field + ME1/PIP, FL1/QUI, CAL1/MOL) with
# PoolSNP: max-cov 0.99, min-cov 10, min-count 20, min-freq 0.05, miss-frac
# 0.2, base-quality 20. Then keeps sites with no missing pool
# (bcftools view -e 'F_MISSING > 0'). Third job of submit_v4.sh.
#
# Step: 4, variant calling (cluster/README.md)
# Manuscript: Methods "Read processing, alignment, and variant calling";
#             Supplementary Methods "Read processing, alignment, variant
#             calling and software versions"
#
# Inputs:  pools27_v4.sanitized.mpileup.gz   (mpileup_v4_merge.sh)
#          bamlist_v4.txt, locations_v4.txt  (pool order; name check only)
#          reference: REF_FASTA_SANITIZED    (02_mpileup/sanitize_headers.sh)
# Outputs: PoolSNP_out_v4/poolsnp.vcf.gz           (raw, <= 20% missing)
#          PoolSNP_out_v4/poolsnp.complete.vcf.gz  (no missing pool) + .tbi
# =============================================================================
#SBATCH --job-name=poolsnp_v4
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH -o logs/poolsnp_v4.%j.out
#SBATCH -e logs/poolsnp_v4.%j.err

set -Eeuo pipefail
trap 'echo "[ERROR] line $LINENO status $?"' ERR

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

BASE="${PROJECT_ROOT}"
REF="${REF_FASTA_SANITIZED}"
MPU="${BASE}/pools27_v4.sanitized.mpileup.gz"
BAMLIST="${BASE}/bamlist_v4.txt"
LOCATIONS="${BASE}/locations_v4.txt"

source "${CONDA_SH}"
conda activate "${ENV_POOLSNP}"

# samtools and bcftools are not in the poolsnp conda env; load from cluster
# modules. Loaded inside the script so neither depends on submit-time state.
module load samtools
module load bcftools

which PoolSNP.sh
which samtools  || { echo "[ERROR] samtools not on PATH after module load";  exit 1; }
which bcftools || { echo "[ERROR] bcftools not on PATH after module load"; exit 1; }
parallel --version

[[ -s "$MPU" && -s "$REF" ]] || { echo "[ERROR] missing input(s)"; exit 1; }

# --- Re-validate bamlist <-> locations alignment (basename only) ------------
# @RG SM check happens at mpileup time; here we just re-confirm the cheap one
# in case the files were edited between mpileup and PoolSNP submission.
i=0
while IFS= read -r bam <&3 && IFS= read -r want <&4; do
  i=$((i+1))
  bn=$(basename "$bam" .filtered.bam)
  [[ "$bn" == "$want" ]] || {
    echo "[ERROR] line $i: BAM basename '$bn' != locations entry '$want'"; exit 1
  }
done 3<"$BAMLIST" 4<"$LOCATIONS"
echo "[INFO] re-validated $i BAM/name pairs"

cd "$BASE"
NAMES=$(paste -sd, locations_v4.txt | tr -d ' \t\r\n')

# Cross-check: names count == mpileup sample count
ns_mpu=$(zcat "$MPU" | awk 'NR==1{v=(NF-3)/3} END{print v}')
ns_names=$(echo "$NAMES" | awk -F, '{print NF}')
echo "[INFO] mpileup samples = $ns_mpu ; names provided = $ns_names"
if [[ "$ns_mpu" -ne "$ns_names" || "$ns_mpu" -ne 27 ]]; then
  echo "[ERROR] sample-count mismatch (expected 27)"; exit 1
fi

# --- PoolSNP parameters -----------------------------------------------------
MINCOV=10
MAXCOV=0.99
MINCOUNT=20
MINFREQ=0.05
MISSFRAC=0.2
JOBS=8
BASEQ=20

OUT="${BASE}/PoolSNP_out_v4/poolsnp"
mkdir -p "$(dirname "$OUT")" logs

echo "[INFO] Running PoolSNP v4..."
bash "$(which PoolSNP.sh)" \
  mpileup="$MPU" \
  reference="$REF" \
  names="$NAMES" \
  max-cov="$MAXCOV" \
  min-cov="$MINCOV" \
  min-count="$MINCOUNT" \
  min-freq="$MINFREQ" \
  miss-frac="$MISSFRAC" \
  jobs="$JOBS" \
  base-quality="$BASEQ" \
  badsites=1 \
  allsites=0 \
  output="$OUT"

VCF_RAW="${OUT}.vcf.gz"
[[ -s "$VCF_RAW" ]] || { echo "[ERROR] expected output not found: $VCF_RAW"; exit 1; }
n_raw=$(zcat "$VCF_RAW" | grep -vc '^#')
echo "[INFO] PoolSNP raw VCF: $VCF_RAW  ($n_raw biallelic SNPs, ≤20% missingness)"

# --- Post-filter: no-missing-data subset ------------------------------------
# bcftools F_MISSING is fraction of samples with missing GT (./.). Requiring
# F_MISSING == 0 keeps only sites with complete AF data across all 27 pools.
VCF_COMPLETE="${OUT}.complete.vcf.gz"
echo "[INFO] Filtering to complete-data subset -> $VCF_COMPLETE"
bcftools view -e 'F_MISSING > 0' "$VCF_RAW" -Oz -o "$VCF_COMPLETE"
bcftools index -t "$VCF_COMPLETE"
n_complete=$(zcat "$VCF_COMPLETE" | grep -vc '^#')
echo "[INFO] complete-data VCF: $VCF_COMPLETE  ($n_complete SNPs)"

echo "[✓] PoolSNP v4 finished."
echo "    raw      : $VCF_RAW           ($n_raw SNPs, ≤20% missing)"
echo "    complete : $VCF_COMPLETE      ($n_complete SNPs, 0 missing)"
