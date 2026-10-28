#!/bin/bash
# =============================================================================
# cluster/06_depth/glmm/03_intersect_cds.slurm
# =============================================================================
# SLURM wrapper for 03_intersect_cds.sh (callable positions restricted to the
# merged CDS of each gene in the lifted annotation).
#
# Step: 8, callable depth (cluster/README.md)
# Inputs:  data/glmm_v4/callable_depth_v4.tsv.gz (02_extract_callable_depth.slurm)
#          reference_genomes/quinq_on_pip.gff3_polished (Liftoff annotation;
#          deposit data/quinq_on_pip.gff3_polished)
# Outputs: data/glmm_v4/cds_per_gene_merged.bed, callable_cds_depth.tsv.gz
# =============================================================================
#SBATCH --job-name=glmm_03_intersect
#SBATCH --output=logs/glmm_03_intersect_%j.out
#SBATCH --error=logs/glmm_03_intersect_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -Eeuo pipefail
cd "$SLURM_SUBMIT_DIR"

module load bedtools

mkdir -p data/glmm_v4 logs

bash scripts/06_depth/glmm/03_intersect_cds.sh \
  data/glmm_v4/callable_depth_v4.tsv.gz \
  reference_genomes/quinq_on_pip.gff3_polished \
  data/glmm_v4
