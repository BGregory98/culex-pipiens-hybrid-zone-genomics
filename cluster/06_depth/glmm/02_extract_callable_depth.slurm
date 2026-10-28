#!/bin/bash
# =============================================================================
# cluster/06_depth/glmm/02_extract_callable_depth.slurm
# =============================================================================
# SLURM wrapper for 02_extract_callable_depth.py (pass 2 over the v4 pileup:
# per-pool depth at every position with depth >= 10 in all 27 pools and no
# pool above its 99th-percentile cutoff).
#
# Step: 8, callable depth (cluster/README.md)
# Inputs:  pools27_v4.sanitized.mpileup.gz, data/metadata/locations_v4.txt,
#          data/glmm_v4/pool_chrom_cutoffs_v4.tsv (01_compute_percentile_cutoffs.slurm)
# Outputs: data/glmm_v4/callable_depth_v4.tsv.gz
# =============================================================================
#SBATCH --job-name=glmm_02_callable
#SBATCH --output=logs/glmm_02_callable_%j.out
#SBATCH --error=logs/glmm_02_callable_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

set -Eeuo pipefail
cd "$SLURM_SUBMIT_DIR"

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

source miniconda3/etc/profile.d/conda.sh
conda activate "${ENV_POOLSNP}"

mkdir -p data/glmm_v4 logs

python scripts/06_depth/glmm/02_extract_callable_depth.py \
  --mpileup   pools27_v4.sanitized.mpileup.gz \
  --locations data/metadata/locations_v4.txt \
  --cutoffs   data/glmm_v4/pool_chrom_cutoffs_v4.tsv \
  --min-cov   10 \
  --out       data/glmm_v4/callable_depth_v4.tsv.gz
