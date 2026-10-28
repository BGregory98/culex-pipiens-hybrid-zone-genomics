#!/bin/bash
# =============================================================================
# cluster/06_depth/glmm/04_aggregate_cds_depth.slurm
# =============================================================================
# SLURM wrapper for 04_aggregate_cds_depth.R with --pool-set field (default):
# per-gene CDS callable depth for the 24 field pools and the genome-wide
# callable-depth offset.
#
# Step: 8, callable depth (cluster/README.md)
# Inputs:  data/glmm_v4/callable_cds_depth.tsv.gz (03_intersect_cds.slurm),
#          data/glmm_v4/callable_depth_v4.tsv.gz, data/metadata/sample_categories.csv
# Outputs: data/glmm_v4/cds_depth_per_gene_pool.tsv,
#          data/glmm_v4/genome_callable_depth_per_pool.tsv
# =============================================================================
#SBATCH --job-name=glmm_04_aggregate
#SBATCH --output=logs/glmm_04_aggregate_%j.out
#SBATCH --error=logs/glmm_04_aggregate_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=1

set -Eeuo pipefail
cd "$SLURM_SUBMIT_DIR"
module load R/4.5.3

mkdir -p data/glmm_v4 logs

Rscript scripts/06_depth/glmm/04_aggregate_cds_depth.R \
  --callable-cds      data/glmm_v4/callable_cds_depth.tsv.gz \
  --callable-all      data/glmm_v4/callable_depth_v4.tsv.gz \
  --sample-categories data/metadata/sample_categories.csv \
  --min-positions     50 \
  --out-dir           data/glmm_v4
