#!/bin/bash
# =============================================================================
# cluster/04_baypass/run_baypass_filter.sh
# =============================================================================
# SLURM wrapper for filter_baypass_results.R: joins the core-run XtXst with
# the IS-run Bayes factors and keeps every SNP x covariate row that passes
# either the POD XtXst 99% threshold or BF > 10 dB. Finishes in 5-15 min.
#
# Usage: sbatch run_baypass_filter.sh <baypass_dir> [<outdir>]
#   <baypass_dir>  directory with the core and IS outputs and the XtX
#                  threshold RDS files (pod_03_apply_threshold.R)
#   <outdir>       directory for the output CSVs [default: <baypass_dir>]
#
# Step: 6, BayPass (cluster/README.md)
# Manuscript: Methods "Genome-environment association with BayPass"
#
# Inputs:  <baypass_dir>/baypass_core_envsubset_summary_pi_xtx.out,
#          baypass_STD_snp_ids.txt, baypass_env_IS_summary_betai_reg.out,
#          XtX_minor_T.rds, XtX_major_T.rds
#          data/metadata/sample_categories.csv,
#          data/metadata/site_coords_with_spatial_data.csv
# Outputs: <outdir>/filtered_baypass_data.csv (read by
#          R/03_baypass/01_baypass_outliers.R), strong_plus_betas.csv
# =============================================================================
#SBATCH --job-name=bp_filter
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --time=0-02:00:00
#SBATCH -o logs/bp_filter.out
#SBATCH -e logs/bp_filter.err

set -euo pipefail

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

BAYPASS_DIR=${1:?Usage: sbatch run_baypass_filter.sh <baypass_dir> [<outdir>]}
OUTDIR=${2:-$BAYPASS_DIR}

METADATA_DIR=$PROJECT_ROOT/data/metadata

cd "$PROJECT_ROOT"
module load R/4.5.3
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.5"

echo "=== Job start: $(date) ==="
echo "Host:        $(hostname)"
echo "Working dir: $(pwd)"
echo "BayPass dir: $BAYPASS_DIR"
echo "Output dir:  $OUTDIR"
echo "R version:   $(Rscript -e 'cat(R.version.string)')"
echo ""

Rscript scripts/04_baypass/filter_baypass_results.R \
  --xtx        "$BAYPASS_DIR/baypass_core_envsubset_summary_pi_xtx.out" \
  --snp_ids    "$BAYPASS_DIR/baypass_STD_snp_ids.txt" \
  --betai_reg  "$BAYPASS_DIR/baypass_env_IS_summary_betai_reg.out" \
  --xtx_minor  "$BAYPASS_DIR/XtX_minor_T.rds" \
  --xtx_major  "$BAYPASS_DIR/XtX_major_T.rds" \
  --categories "$METADATA_DIR/sample_categories.csv" \
  --env        "$METADATA_DIR/site_coords_with_spatial_data.csv" \
  --outdir     "$OUTDIR"

echo ""
echo "=== Job end: $(date) ==="
ls -lh "$OUTDIR/filtered_baypass_data.csv" "$OUTDIR/strong_plus_betas.csv" 2>/dev/null || true
