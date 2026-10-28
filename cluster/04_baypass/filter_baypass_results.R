#!/usr/bin/env Rscript
# =============================================================================
# cluster/04_baypass/filter_baypass_results.R
# =============================================================================
# Joins the core-run XtXst (one value per SNP) with the IS-run Bayes factors
# (one per SNP x covariate), keeps rows with XtXst above the POD 99% threshold
# or BF above 10 dB, and labels each row's XtXst and BF as minor, major or
# none. BayPass reports BF in deciban units, BF.dB. = 10 log10(BF), so the
# 10 dB and 30 dB cut-offs are BF = 10 and BF = 1,000. Submitted through
# run_baypass_filter.sh.
#
# Step: 6, BayPass (cluster/README.md)
# Manuscript: Methods "Genome-environment association with BayPass"
#
# Inputs:  --xtx, --snp_ids, --betai_reg (baypass_core_envsubset.sh,
#          prepare_baypass_input_v4.R, baypass_env_IS.sh); --xtx_minor,
#          --xtx_major (pod_03_apply_threshold.R); --categories, --env
#          (data/metadata/)
# Outputs: see below; filtered_baypass_data.csv is read by
#          R/03_baypass/01_baypass_outliers.R (deposit
#          results/baypass/v4/filtered_baypass_data.csv)
#
# Usage:
#   Rscript filter_baypass_results.R \
#     --xtx        <core_outprefix>_summary_pi_xtx.out \
#     --snp_ids    baypass_STD_snp_ids.txt \
#     --betai_reg  <is_outprefix>_summary_betai_reg.out \
#     --xtx_minor  XtX_minor_T.rds \
#     --xtx_major  XtX_major_T.rds \
#     --categories <metadata_dir>/sample_categories.csv \
#     --env        <metadata_dir>/site_coords_with_spatial_data.csv \
#     --outdir     <directory for filtered outputs>
#
# Outputs (to --outdir):
#   strong_plus_betas.csv      all SNP x predictor rows with BF > 10 dB
#   filtered_baypass_data.csv  rows passing XtXst > XtX_minor_T or
#                              BF.dB. > 10, with XtX_outlier_type and
#                              BF_outlier_type labels
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(optparse)
})

# ---- CLI args -------------------------------------------------------------

opt_list <- list(
  make_option("--xtx",        type = "character", default = NULL,
              help = "Path to core run's _summary_pi_xtx.out file"),
  make_option("--snp_ids",    type = "character", default = NULL,
              help = "Path to baypass_STD_snp_ids.txt (SNP IDs aligned with XtX)"),
  make_option("--betai_reg",  type = "character", default = NULL,
              help = "Path to IS run's _summary_betai_reg.out file"),
  make_option("--xtx_minor",  type = "character", default = NULL,
              help = "Path to XtX_minor_T.rds (99% POD-calibrated threshold)"),
  make_option("--xtx_major",  type = "character", default = NULL,
              help = "Path to XtX_major_T.rds (99.9% POD-calibrated threshold)"),
  make_option("--categories", type = "character", default = NULL,
              help = "Path to sample_categories.csv"),
  make_option("--env",        type = "character", default = NULL,
              help = "Path to site_coords_with_spatial_data.csv"),
  make_option("--outdir",     type = "character", default = NULL,
              help = "Directory for output CSVs")
)
opt <- parse_args(OptionParser(option_list = opt_list))

required <- c("xtx", "snp_ids", "betai_reg", "xtx_minor", "xtx_major",
              "categories", "env", "outdir")
missing  <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing) > 0) {
  stop("Missing required args: --", paste(missing, collapse = ", --"),
       "\nRun with --help for usage.")
}

XTX_PATH       <- opt$xtx
SNP_IDS_PATH   <- opt$snp_ids
BETAI_REG_PATH <- opt$betai_reg
XTX_MINOR_PATH <- opt$xtx_minor
XTX_MAJOR_PATH <- opt$xtx_major
CATEGORIES     <- opt$categories
ENV_PATH       <- opt$env
OUTDIR         <- opt$outdir

# Bayes-factor thresholds in deciban units (dB = 10 log10 BF), fixed rather
# than calibrated on the data. 10 dB = BF 10 ("strong" on Jeffreys' scale);
# 30 dB = BF 1,000 (in Jeffreys' "decisive" range, which starts at 20 dB).
# BF_MAJOR is used only for the BF_outlier_type label, not for the row filter.
BF_MINOR <- 10    # dB (BF = 10)
BF_MAJOR <- 30    # dB (BF = 1,000)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---- Sanity check: required files exist ----------------------------------

required_files <- c(XTX_PATH, SNP_IDS_PATH, BETAI_REG_PATH,
                    XTX_MINOR_PATH, XTX_MAJOR_PATH, CATEGORIES, ENV_PATH)
missing_f <- required_files[!file.exists(required_files)]
if (length(missing_f) > 0) {
  stop("Missing required input file(s):\n  ",
       paste(missing_f, collapse = "\n  "))
}

# ---- 1. Sample / env metadata --------------------------------------------

sample_info <- read.csv(CATEGORIES)
field_names <- sort(sample_info$sample[sample_info$sample_type == "field"])

env_raw  <- read.csv(ENV_PATH)
env_raw  <- env_raw %>% dplyr::select(-c(lat, lon)) %>% rename(sample = site)
env_data <- env_raw %>%
  filter(sample %in% field_names) %>%
  arrange(match(sample, field_names))

cat(sprintf("Field samples: %d\n", length(field_names)))

# ---- 2. XtX statistics with SNP IDs --------------------------------------

xtx <- data.table::fread(XTX_PATH, header = TRUE, data.table = FALSE, check.names = TRUE)
baypass_STD_snp_ids <- read.table(SNP_IDS_PATH)
xtx <- cbind(baypass_STD_snp_ids, xtx) %>%
  separate(
    V1,
    into    = c("chromosome", "position"),
    sep     = "_(?=[0-9eE.+-]+$)",   # underscore before a numeric/scientific block
    remove  = TRUE,
    convert = TRUE
  )

cat(sprintf("Loaded XtX table: %d SNPs\n", nrow(xtx)))

# ---- 3. Beta coefficients (BFs per SNP × predictor) ----------------------

betas <- data.table::fread(BETAI_REG_PATH, header = TRUE, data.table = FALSE, check.names = TRUE)
betas$COVARIABLE <- as.factor(betas$COVARIABLE)

# Map BayPass numeric COVARIABLE indices back to predictor names from env_data
covariable_names      <- names(env_data[, -1])
covariable_dictionary <- data.frame(
  predictor  = as.factor(covariable_names),
  COVARIABLE = as.factor(seq_along(covariable_names))
)

betas        <- full_join(betas, covariable_dictionary, by = "COVARIABLE")
betas$BF_raw <- 10^(betas$BF.dB. / 10)

# Save the small "strong_plus_betas" table. This compares the raw BF with the
# dB threshold, which is exact only because 10 dB = BF 10: the row set is
# BF.dB. > 10, the same as the filter below.
strong_plus_betas <- betas %>% filter(BF_raw > BF_MINOR)
write.csv(
  strong_plus_betas,
  file.path(OUTDIR, "strong_plus_betas.csv"),
  row.names = FALSE
)
cat(sprintf("Wrote strong_plus_betas.csv (%d rows)\n", nrow(strong_plus_betas)))

# ---- 4. Load POD-calibrated XtX thresholds -------------------------------

XtX_minor_T <- readRDS(XTX_MINOR_PATH)
XtX_major_T <- readRDS(XTX_MAJOR_PATH)

cat("\nThresholds:\n")
cat(sprintf("  XtX minor (99%%):   %.4f\n", as.numeric(XtX_minor_T)))
cat(sprintf("  XtX major (99.9%%): %.4f\n", as.numeric(XtX_major_T)))
cat(sprintf("  BF  minor (dB):    %.4f\n", BF_MINOR))
cat(sprintf("  BF  major (dB):    %.4f\n", BF_MAJOR))

# ---- 5. Filter-first-then-join (memory-efficient) ------------------------
#
# A full_join(betas, xtx) before filtering builds a ~43M-row, multi-GB
# intermediate. Instead, identify MRKs surviving the filter BEFORE joining:
#   - MRKs where XtXst > XtX_minor_T  (XtX is per-MRK, so any row for that
#     MRK passes the OR-clause regardless of BF)
#   - MRKs where any (MRK, predictor) row has BF.dB. > BF_MINOR
# Subset both tables to just those MRKs, join the small subsets, then apply
# the original per-row filter.
#
# Output is row-for-row identical to:
#   full_join(betas, xtx) %>% filter(XtXst > XtX_minor_T | BF.dB. > BF_MINOR)

t_filter_start <- Sys.time()

mrks_pass_xtx <- xtx %>% filter(XtXst > XtX_minor_T) %>% pull(MRK)
mrks_pass_bf  <- strong_plus_betas$MRK   # already BF.dB. > BF_MINOR
keep_mrks     <- union(unique(mrks_pass_xtx), unique(mrks_pass_bf))

cat(sprintf("\nMRKs passing XtX threshold: %d\n", length(unique(mrks_pass_xtx))))
cat(sprintf("MRKs passing BF threshold:  %d\n", length(unique(mrks_pass_bf))))
cat(sprintf("Union (MRKs to keep):       %d\n", length(keep_mrks)))

xtx_sub   <- xtx   %>% filter(MRK %in% keep_mrks)
betas_sub <- betas %>% filter(MRK %in% keep_mrks)

cat(sprintf("xtx_sub rows:   %d\nbetas_sub rows: %d\n",
            nrow(xtx_sub), nrow(betas_sub)))

# Join is now small
bp_joined <- full_join(betas_sub, xtx_sub, by = "MRK") %>%
  dplyr::select(
    MRK, chromosome, position,
    M_P, SD_P, M_XtX, SD_XtX, XtXst, log10.1.pval.,
    COVARIABLE, predictor,
    M_Pearson, SD_Pearson, Beta_is, SD_Beta_is, eBPis, BF.dB., BF_raw
  )

# Apply the per-row filter
bp_outliers_raw <- bp_joined %>%
  filter(XtXst > XtX_minor_T | BF.dB. > BF_MINOR)

bp_outliers_raw <- bp_outliers_raw %>%
  mutate(
    XtX_outlier_type = case_when(
      XtXst > XtX_major_T                        ~ 'major',
      XtXst > XtX_minor_T & XtXst <= XtX_major_T ~ 'minor',
      TRUE                                       ~ 'none'
    ),
    BF_outlier_type = case_when(
      BF.dB. > BF_MAJOR                         ~ 'major',
      BF.dB. > BF_MINOR & BF.dB. <= BF_MAJOR    ~ 'minor',
      TRUE                                      ~ 'none'
    )
  )

t_filter_end <- Sys.time()
cat(sprintf("\nFilter+join elapsed: %.1f seconds\n",
            as.numeric(difftime(t_filter_end, t_filter_start, units = "secs"))))

# ---- 6. Save filtered output ---------------------------------------------

out_path <- file.path(OUTDIR, "filtered_baypass_data.csv")
write.csv(bp_outliers_raw, out_path, row.names = FALSE)

cat(sprintf(
  "\nWrote %s\n  rows: %d\n  unique SNPs: %d\n",
  out_path,
  nrow(bp_outliers_raw),
  length(unique(bp_outliers_raw$MRK))
))

cat("\nDone. Outputs for the R stage:\n")
cat(sprintf("  %s\n", file.path(OUTDIR, "filtered_baypass_data.csv")))
cat(sprintf("  %s\n", file.path(OUTDIR, "strong_plus_betas.csv")))
cat("Next: R/03_baypass/01_baypass_outliers.R.\n")
