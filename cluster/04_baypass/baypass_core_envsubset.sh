#!/bin/bash
# =============================================================================
# cluster/04_baypass/baypass_core_envsubset.sh
# =============================================================================
# BayPass 3.1 core model on the 24 field pools (Pool-seq mode, haploid pool
# size 80): estimates the population covariance matrix (Omega), the
# (a_pi, b_pi) parameters and calibrated XtXst. -nval 1000 -thin 25
# -burnin 5000; pilot runs left at the defaults (10 x 100 iterations, as the
# run log records). Ran ~8 days on 32 threads.
#
# Step: 6, BayPass (cluster/README.md)
# Manuscript: Methods "Genome-environment association with BayPass"
#
# Inputs:  $INDIR/baypass_STD_counts.txt, baypass_poolsize_envsubset.txt
#          (prepare_baypass_input_v4.R)
# Outputs: $INDIR/baypass_core_envsubset_mat_omega.out,
#          _summary_beta_params.out, _summary_pi_xtx.out (and BayPass's other
#          _summary_* files)
# =============================================================================
#SBATCH --job-name=baypass_core_env
#SBATCH -c 32
#SBATCH --mem=256G
#SBATCH --time=10-00:00:00
#SBATCH -o logs/baypass_core_env.out
#SBATCH -e logs/baypass_core_env.err

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

# Input directory, relative to the project root; override at submission, e.g.
#   INDIR=SequencingData/data/baypass_vX sbatch scripts/04_baypass/baypass_core_envsubset.sh
INDIR=${INDIR:-SequencingData/data/baypass_v4}

cd "${PROJECT_ROOT}"

"$BAYPASS" \
  -gfile        "$INDIR/baypass_STD_counts.txt" \
  -poolsizefile "$INDIR/baypass_poolsize_envsubset.txt" \
  -outprefix    "$INDIR/baypass_core_envsubset" \
  -nval 1000 \
  -thin 25 \
  -burnin 5000 \
  -nthreads "${SLURM_CPUS_PER_TASK:-8}"
