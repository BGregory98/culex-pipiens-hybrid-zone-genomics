#!/usr/bin/env Rscript
# =============================================================================
# cluster/04_baypass/pod_03_apply_threshold.R
# =============================================================================
# Takes the 99%, 99.9% and 99.99% quantiles of the POD XtXst distribution as
# the minor, major and extreme outlier thresholds (published values 71.1 and
# 106.6 for minor and major) and, given the real-data XtXst, counts the SNPs
# above each.
#
# Step: 6, BayPass (cluster/README.md)
# Manuscript: Methods "Genome-environment association with BayPass"
#
# Inputs:  --pod_xtx (pod_02_run_baypass.sh); --real_xtx
#          (baypass_core_envsubset.sh)
# Outputs: <outdir>/XtX_minor_T.rds, XtX_major_T.rds (read by
#          filter_baypass_results.R), XtX_extreme_T.rds,
#          pod_threshold_summary.csv
#
# Usage:
#   Rscript pod_03_apply_threshold.R \
#     --pod_xtx  <pod_outprefix>_summary_pi_xtx.out \
#     --outdir   <directory for threshold RDS files> \
#     [--real_xtx <real_run_outprefix>_summary_pi_xtx.out] \
#     [--empirical_thresholds_dir <dir with legacy XtX_99_T.rds etc.>]
#
# Published run (from the project root, R/4.5.3 module):
#   Rscript scripts/04_baypass/pod_03_apply_threshold.R \
#     --pod_xtx  SequencingData/data/baypass_v4/pod/baypass_pod_summary_pi_xtx.out \
#     --real_xtx SequencingData/data/baypass_v4/baypass_core_envsubset_summary_pi_xtx.out \
#     --outdir   SequencingData/data/baypass_v4
#
# Outputs (to --outdir):
#   XtX_minor_T.rds            99%    POD-calibrated XtXst threshold
#   XtX_major_T.rds            99.9%  POD-calibrated XtXst threshold
#   XtX_extreme_T.rds          99.99% POD-calibrated XtXst threshold (not used downstream)
#   pod_threshold_summary.csv  threshold values + outlier counts (if --real_xtx given)
#
# --empirical_thresholds_dir prints a comparison with empirical (genome-wide)
# quantile thresholds saved as XtX_{99,999,9999}_T.rds; not used for the paper.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# ---- CLI args -------------------------------------------------------------

opt_list <- list(
  make_option("--pod_xtx",                  type = "character", default = NULL,
              help = "Path to POD run's _summary_pi_xtx.out file"),
  make_option("--outdir",                   type = "character", default = NULL,
              help = "Directory for output threshold RDS files + summary CSV"),
  make_option("--real_xtx",                 type = "character", default = NULL,
              help = "Optional: real-data _summary_pi_xtx.out for outlier counts at POD thresholds"),
  make_option("--empirical_thresholds_dir", type = "character", default = NULL,
              help = "Optional: dir with legacy XtX_{99,999,9999}_T.rds for side-by-side comparison")
)
opt <- parse_args(OptionParser(option_list = opt_list))

required   <- c("pod_xtx", "outdir")
required_r <- gsub("-", "_", required)
missing    <- required[vapply(required_r, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing) > 0) {
  stop("Missing required args: --", paste(missing, collapse = ", --"),
       "\nRun with --help for usage.")
}

POD_XTX_PATH <- opt$pod_xtx
OUTDIR       <- opt$outdir
REAL_XTX     <- opt$real_xtx
EMP_DIR      <- opt$empirical_thresholds_dir

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---- Read POD XtX ---------------------------------------------------------

if (!file.exists(POD_XTX_PATH)) {
  stop("POD XtX file not found: ", POD_XTX_PATH,
       "\nDid pod_02_run_baypass.sh complete?")
}
pod_xtx <- fread(POD_XTX_PATH)
cat("POD SNPs:", nrow(pod_xtx), "\n")
cat("POD XtXst summary:\n"); print(summary(pod_xtx$XtXst))

# ---- Compute POD-calibrated thresholds -----------------------------------

XtX_minor_T   <- quantile(pod_xtx$XtXst, 0.99,   na.rm = TRUE)
XtX_major_T   <- quantile(pod_xtx$XtXst, 0.999,  na.rm = TRUE)
XtX_extreme_T <- quantile(pod_xtx$XtXst, 0.9999, na.rm = TRUE)

cat(sprintf("\nPOD-calibrated XtX thresholds:\n"))
cat(sprintf("  99%%    (minor)  : %.4f\n", XtX_minor_T))
cat(sprintf("  99.9%%  (major)  : %.4f\n", XtX_major_T))
cat(sprintf("  99.99%% (extreme): %.4f\n", XtX_extreme_T))

# ---- Save (with v4-pipeline-doc naming) ----------------------------------

saveRDS(XtX_minor_T,   file.path(OUTDIR, "XtX_minor_T.rds"))
saveRDS(XtX_major_T,   file.path(OUTDIR, "XtX_major_T.rds"))
saveRDS(XtX_extreme_T, file.path(OUTDIR, "XtX_extreme_T.rds"))
cat("\nWrote: XtX_minor_T.rds, XtX_major_T.rds, XtX_extreme_T.rds to", OUTDIR, "\n")

# ---- Outlier counts in real data (optional) -------------------------------

n_real_minor   <- NA_integer_
n_real_major   <- NA_integer_
n_real_extreme <- NA_integer_

if (!is.null(REAL_XTX)) {
  if (!file.exists(REAL_XTX)) {
    warning("--real_xtx given but file not found: ", REAL_XTX,
            "\nSkipping outlier counts.")
  } else {
    real_xtx <- fread(REAL_XTX)
    n_total  <- nrow(real_xtx)
    n_real_minor   <- sum(real_xtx$XtXst > XtX_minor_T,   na.rm = TRUE)
    n_real_major   <- sum(real_xtx$XtXst > XtX_major_T,   na.rm = TRUE)
    n_real_extreme <- sum(real_xtx$XtXst > XtX_extreme_T, na.rm = TRUE)
    cat(sprintf("\n--- Outlier counts in real data (n = %d SNPs) ---\n", n_total))
    cat(sprintf("  > minor   (99%%)   : %d (%.3f%%)\n", n_real_minor,   100 * n_real_minor   / n_total))
    cat(sprintf("  > major   (99.9%%) : %d (%.4f%%)\n", n_real_major,   100 * n_real_major   / n_total))
    cat(sprintf("  > extreme (99.99%%): %d (%.5f%%)\n", n_real_extreme, 100 * n_real_extreme / n_total))
  }
}

# ---- Legacy empirical-threshold comparison (optional) ---------------------

if (!is.null(EMP_DIR)) {
  emp_paths <- list(
    minor   = file.path(EMP_DIR, "XtX_99_T.rds"),
    major   = file.path(EMP_DIR, "XtX_999_T.rds"),
    extreme = file.path(EMP_DIR, "XtX_9999_T.rds")
  )
  if (all(sapply(emp_paths, file.exists))) {
    emp <- lapply(emp_paths, readRDS)
    cat("\n--- Empirical vs POD-calibrated comparison ---\n")
    cat(sprintf("  Quantile      Empirical       POD-calibrated     Δ\n"))
    cat(sprintf("  99%%        %12.4f       %12.4f   %+8.4f\n",
                emp$minor,   XtX_minor_T,   XtX_minor_T   - emp$minor))
    cat(sprintf("  99.9%%      %12.4f       %12.4f   %+8.4f\n",
                emp$major,   XtX_major_T,   XtX_major_T   - emp$major))
    cat(sprintf("  99.99%%     %12.4f       %12.4f   %+8.4f\n",
                emp$extreme, XtX_extreme_T, XtX_extreme_T - emp$extreme))
  } else {
    warning("--empirical_thresholds_dir given but expected files not all found in ",
            EMP_DIR, "; skipping comparison.")
  }
}

# ---- Save summary CSV -----------------------------------------------------

summary_dt <- data.table(
  quantile   = c("99%", "99.9%", "99.99%"),
  label      = c("minor", "major", "extreme"),
  threshold  = c(XtX_minor_T, XtX_major_T, XtX_extreme_T),
  n_outliers = c(n_real_minor, n_real_major, n_real_extreme)
)
fwrite(summary_dt, file.path(OUTDIR, "pod_threshold_summary.csv"))
cat("\nWrote: pod_threshold_summary.csv\n")

cat("\nDone.\n")
