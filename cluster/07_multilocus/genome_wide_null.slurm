#!/bin/bash
# =============================================================================
# cluster/07_multilocus/genome_wide_null.slurm
# =============================================================================
# SLURM wrapper for 01_genome_wide_null.R: 10,000 random matched-width
# regions per locus, seed 42 (N_REGIONS and SEED can be overridden from the
# environment). Submitted from the project root:
#   sbatch scripts/07_multilocus/genome_wide_null.slurm
# Ran ~5 h.
#
# Step: 9, genome-wide null (cluster/README.md)
# Inputs/Outputs: see 01_genome_wide_null.R; writes genome_wide_null.rds in
#                 the project root
# =============================================================================
#SBATCH --job-name=ml_null
#SBATCH --output=logs/ml_null_%j.out
#SBATCH --error=logs/ml_null_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

set -Eeuo pipefail

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs SequencingData/results/multilocus

module purge
module load R/4.5.3

# Override from environment (e.g. `N_REGIONS=5000 sbatch ...`)
N_REGIONS=${N_REGIONS:-10000}
SEED=${SEED:-42}

echo "=== genome_wide_null ==="
echo "  n_regions: $N_REGIONS"
echo "  seed:      $SEED"
echo "  Job:       $SLURM_JOB_ID"
echo "  Host:      $(hostname)"
echo "  Started:   $(date)"

Rscript scripts/07_multilocus/01_genome_wide_null.R \
  --base_dir       "$SLURM_SUBMIT_DIR" \
  --vcf            "PoolSNP_out_v4/poolsnp.biallelic.nomiss.vcf.gz" \
  --env_file       "data/metadata/site_coords_with_spatial_data.csv" \
  --out_rds        "genome_wide_null.rds" \
  --n_regions      "$N_REGIONS" \
  --seed           "$SEED"

echo "  Finished:  $(date)"

