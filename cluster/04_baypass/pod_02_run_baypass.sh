#!/bin/bash
# =============================================================================
# cluster/04_baypass/pod_02_run_baypass.sh
# =============================================================================
# BayPass 3.1 core model on the pseudo-observed dataset (POD) simulated by
# pod_01_generate.R, with the real run's MCMC settings and 10 pilot runs of
# 100 iterations. Its XtXst distribution is the neutral null that
# pod_03_apply_threshold.R turns into the outlier thresholds. Ran ~3 h.
#
# Usage: sbatch pod_02_run_baypass.sh <pod_dir> [<out_prefix>]
#   <pod_dir>     directory with G.btapods and poolsize.btapods
#   <out_prefix>  output file prefix [default: baypass_pod]
#
# Step: 6, BayPass (cluster/README.md)
# Manuscript: Methods "Genome-environment association with BayPass"
#
# Inputs:  <pod_dir>/G.btapods, poolsize.btapods (pod_01_generate.R)
# Outputs: <pod_dir>/<out_prefix>_summary_pi_xtx.out (and BayPass's other
#          output files)
# =============================================================================
#SBATCH --job-name=baypass_pod
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH -o logs/baypass_pod.out
#SBATCH -e logs/baypass_pod.err

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

PODDIR=${1:?Usage: sbatch pod_02_run_baypass.sh <pod_dir> [<out_prefix>]}
OUT_PREFIX=${2:-baypass_pod}

cd "$PODDIR" || { echo "ERROR: cannot cd to $PODDIR"; exit 1; }

# Sanity checks on inputs from pod_01_generate.R
if [[ ! -f G.btapods ]]; then
    echo "ERROR: G.btapods not found in $PODDIR. Run pod_01_generate.R first." >&2
    exit 1
fi
if [[ ! -f poolsize.btapods ]]; then
    echo "ERROR: poolsize.btapods not found in $PODDIR. Run pod_01_generate.R first." >&2
    exit 1
fi

echo "Starting BayPass POD run at $(date)"
echo "Input dir:     $PODDIR"
echo "G.btapods:     $(wc -l < G.btapods) lines"
echo "Output prefix: $OUT_PREFIX"

"$BAYPASS" \
    -gfile        G.btapods \
    -poolsizefile poolsize.btapods \
    -nval 1000\
    -thin 25\
    -burnin 5000\
    -npilot 10\
    -pilotlength 100\
    -outprefix    "$OUT_PREFIX" \
    -nthreads     "${SLURM_CPUS_PER_TASK:-8}"

echo "Finished at $(date)"
