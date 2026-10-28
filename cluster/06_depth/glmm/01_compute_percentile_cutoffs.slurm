#!/bin/bash
# =============================================================================
# cluster/06_depth/glmm/01_compute_percentile_cutoffs.slurm
# =============================================================================
# SLURM wrapper for 01_compute_percentile_cutoffs.py (pass 1 over the v4
# pileup: per-pool, per-chromosome 99th-percentile depth). Submitted from the
# project root, first of the four-job chain (see README.md here).
#
# Step: 8, callable depth (cluster/README.md)
# Inputs:  pools27_v4.sanitized.mpileup.gz (_toplevel/mpileup_v4_merge.sh)
#          data/metadata/locations_v4.txt
# Outputs: data/glmm_v4/pool_chrom_cutoffs_v4.tsv
# =============================================================================
#SBATCH --job-name=glmm_01_cutoffs
#SBATCH --output=logs/glmm_01_cutoffs_%j.out
#SBATCH --error=logs/glmm_01_cutoffs_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

set -Eeuo pipefail
cd "$SLURM_SUBMIT_DIR"

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

source miniconda3/etc/profile.d/conda.sh
conda activate "${ENV_POOLSNP}"

mkdir -p data/glmm_v4 logs

python scripts/06_depth/glmm/01_compute_percentile_cutoffs.py \
  --mpileup   pools27_v4.sanitized.mpileup.gz \
  --locations data/metadata/locations_v4.txt \
  --out       data/glmm_v4/pool_chrom_cutoffs_v4.tsv
