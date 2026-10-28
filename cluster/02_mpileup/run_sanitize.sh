#!/bin/bash
# =============================================================================
# cluster/02_mpileup/run_sanitize.sh
# =============================================================================
# SLURM wrapper for sanitize_headers.sh, submitted from the directory that
# holds both scripts.
#
# Step: 4, variant calling (reference preparation; cluster/README.md)
# Inputs/Outputs: see sanitize_headers.sh
# =============================================================================
#SBATCH --job-name=sanitize
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH -o logs/sanitize.%j.out
#SBATCH -e logs/sanitize.%j.err

bash sanitize_headers.sh

