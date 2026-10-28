#!/bin/bash
# =============================================================================
# cluster/05_windowed_stats/prepare_pooldata.slurm
# =============================================================================
# SLURM wrapper for prepare_pooldata.R. Submitted from the project root:
#   sbatch --export=WIN_SIZE=1000000 scripts/05_windowed_stats/prepare_pooldata.slurm
#
# Step: 7, windowed f3 (cluster/README.md)
# Inputs/Outputs: see prepare_pooldata.R
# =============================================================================
#SBATCH --job-name=prep_pooldata
#SBATCH --output=logs/prep_pooldata_%j.out
#SBATCH --error=logs/prep_pooldata_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
cd "$SLURM_SUBMIT_DIR"
module load R/4.5.3
Rscript scripts/05_windowed_stats/prepare_pooldata.R ${WIN_SIZE:-1000000}
