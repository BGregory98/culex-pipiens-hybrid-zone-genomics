#!/bin/bash
# =============================================================================
# cluster/04_baypass/baypass_env_IS.sh
# =============================================================================
# BayPass 3.1 standard covariate model in importance-sampling mode, with Omega
# fixed from the core run and the eight scaled covariates (elevation, pis,
# tree_canopy, pd, wtd, bio05, bio07, bio18). Gives the per-SNP x covariate
# Bayes factors (BF in dB) and beta coefficients. Same MCMC settings as the
# core run.
#
# Step: 6, BayPass (cluster/README.md)
# Manuscript: Methods "Genome-environment association with BayPass"
#
# Inputs:  $INDIR/baypass_STD_counts.txt, baypass_poolsize_envsubset.txt,
#          baypass_env_matrix.txt (prepare_baypass_input_v4.R)
#          $INDIR/baypass_core_envsubset_mat_omega.out (baypass_core_envsubset.sh)
# Outputs: $INDIR/baypass_env_IS_summary_betai_reg.out (and BayPass's other
#          _summary_* files)
# =============================================================================
#SBATCH --job-name=baypass_env_IS
#SBATCH -c 32
#SBATCH --mem=256G
#SBATCH --time=10-00:00:00
#SBATCH -o logs/baypass_env_IS.out
#SBATCH -e logs/baypass_env_IS.err

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

# Input directory, relative to the project root; override at submission, e.g.
#   INDIR=SequencingData/data/baypass_vX sbatch scripts/04_baypass/baypass_env_IS.sh
INDIR=${INDIR:-SequencingData/data/baypass_v4}

cd "${PROJECT_ROOT}"

"$BAYPASS" \
  -gfile        "$INDIR/baypass_STD_counts.txt" \
  -poolsizefile "$INDIR/baypass_poolsize_envsubset.txt" \
  -efile        "$INDIR/baypass_env_matrix.txt" \
  -omegafile    "$INDIR/baypass_core_envsubset_mat_omega.out" \
  -outprefix    "$INDIR/baypass_env_IS" \
  -nval 1000 \
  -thin 25 \
  -burnin 5000 \
  -nthreads "${SLURM_CPUS_PER_TASK:-8}"
