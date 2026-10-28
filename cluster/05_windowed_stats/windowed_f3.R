# =============================================================================
# cluster/05_windowed_stats/windowed_f3.R
# =============================================================================
# Computes windowed f3*(focal_pop; ref1, ref2) at a given window size for one
# focal population and reference pairing per invocation (one array task of
# run_windowed_f3.sh). Block-jackknife standard errors use blocks of ~10% of
# each window's SNPs (nsnp.per.bjack.block = n_snps / 10).
#
# Step: 7, windowed f3 (cluster/README.md)
# Manuscript: Methods "Windowed f3 admixture scan"
#
# f3* is the bias-corrected three-population test from poolfstat. Negative
# f3* with significantly negative z indicates focal_pop has admixed ancestry
# from populations related to ref1 and ref2.
#
# Usage:
#   Rscript windowed_f3.R <focal_pop> <ref1> <ref2> <win_size_bp> <out_csv>
#
# Inputs:
#   SequencingData/data/pooldata_v4.rds                (from prepare_pooldata.R)
#   SequencingData/data/windows_tbl_v4_<size>.rds      (from prepare_pooldata.R)
#
# Output:
#   <out_csv>  one row per window: estimate, bjack_mean, bjack_se, z_score,
#              window_chr, window_start, window_end, n_snps, popA (focal),
#              popB (ref1), popC (ref2), win_size_bp
#
# P values and the Benjamini-Hochberg correction are computed after all
# tasks finish, in combine_windowed_outputs.R.
# =============================================================================

suppressPackageStartupMessages({
  library(poolfstat)
  library(data.table)
  library(dplyr)
  library(tibble)
})

# ---- Parse arguments ------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript windowed_f3.R <focal_pop> <ref1> <ref2> <win_size_bp> <out_csv>")
}

focal_pop   <- args[1]
ref1        <- args[2]
ref2        <- args[3]
win_size_bp <- as.integer(args[4])
out_csv     <- args[5]

if (is.na(win_size_bp) || win_size_bp < 1000) {
  stop("win_size_bp must be a positive integer in bp.")
}

triplet <- c(focal_pop, ref1, ref2)
if (length(unique(triplet)) != 3) {
  stop("focal/ref1/ref2 must all differ; got: ",
       paste(triplet, collapse = ", "))
}

format_win_size <- function(n) {
  if (n %% 1e6 == 0) return(paste0(n %/% 1e6, "M"))
  if (n %% 1e3 == 0) return(paste0(n %/% 1e3, "K"))
  return(as.character(n))
}
WIN_SIZE_LABEL <- format_win_size(win_size_bp)

POOLDATA_RDS <- "SequencingData/data/pooldata_v4.rds"
WINDOWS_RDS  <- sprintf("SequencingData/data/windows_tbl_v4_%s.rds", WIN_SIZE_LABEL)

# ---- Sanity checks --------------------------------------------------------

if (!file.exists(POOLDATA_RDS)) {
  stop("pooldata RDS not found at ", POOLDATA_RDS,
       ". Run prepare_pooldata.R first.")
}
if (!file.exists(WINDOWS_RDS)) {
  stop("windows_tbl RDS not found at ", WINDOWS_RDS,
       ". Run: Rscript prepare_pooldata.R ", win_size_bp)
}

# ---- Load preprocessed data -----------------------------------------------

cat("[", format(Sys.time()), "] Loading pooldata + windows (",
    WIN_SIZE_LABEL, ")...\n")
pooldata    <- readRDS(POOLDATA_RDS)
windows_tbl <- readRDS(WINDOWS_RDS)

pool_names <- pooldata@poolnames
for (p in triplet) {
  if (!(p %in% pool_names)) {
    stop("Pop not in pooldata: ", p,
         "\n  Available: ", paste(pool_names, collapse = ", "))
  }
}

p_idx_focal <- which(pool_names == focal_pop)
p_idx_ref1  <- which(pool_names == ref1)
p_idx_ref2  <- which(pool_names == ref2)

cat("  Focal:  ", focal_pop, "(idx", p_idx_focal, ")\n")
cat("  Ref1:   ", ref1, "(idx", p_idx_ref1, ")\n")
cat("  Ref2:   ", ref2, "(idx", p_idx_ref2, ")\n")
cat("  Windows:", nrow(windows_tbl), "\n")

# ---- Subset pools ONCE outside the window loop ---------------------------

cat("[", format(Sys.time()), "] Subsetting pooldata to focal+ref1+ref2...\n")
pooldata_triplet <- pooldata.subset(
  pooldata,
  pool.index = c(p_idx_focal, p_idx_ref1, p_idx_ref2),
  verbose    = FALSE
)
rm(pooldata); gc()

# ---- Loop over windows; collect results in a pre-allocated list ----------

cat("[", format(Sys.time()), "] Running f3* on", nrow(windows_tbl),
    "windows...\n")
results <- vector("list", nrow(windows_tbl))
start_time <- Sys.time()

for (w in seq_len(nrow(windows_tbl))) {
  snp_idx_vec <- windows_tbl$snp_idx[[w]]
  if (length(snp_idx_vec) == 0L) next

  pooldata_window <- pooldata.subset(
    pooldata_triplet,
    snp.index = snp_idx_vec,
    verbose   = FALSE
  )

  fstats_window <- compute.fstats(
    pooldata_window,
    nsnp.per.bjack.block = round(windows_tbl$n_snps[w] / 10),
    computeQmat          = FALSE,
    computeF4            = FALSE,
    verbose              = FALSE
  )
  f3star <- head(as.data.frame(fstats_window@f3star.values), 1)
  names(f3star) <- c("estimate", "bjack_mean", "bjack_se", "z_score")

  results[[w]] <- data.table(
    estimate     = f3star$estimate,
    bjack_mean   = f3star$bjack_mean,
    bjack_se     = f3star$bjack_se,
    z_score      = f3star$z_score,
    window_chr   = as.character(windows_tbl$chromosome[w]),
    window_start = windows_tbl$start[w],
    window_end   = windows_tbl$end[w],
    n_snps       = windows_tbl$n_snps[w],
    popA         = focal_pop,
    popB         = ref1,
    popC         = ref2,
    win_size_bp  = win_size_bp
  )

  if (w %% 50 == 0 || w == nrow(windows_tbl)) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    rate    <- w / elapsed
    eta     <- (nrow(windows_tbl) - w) / rate
    cat(sprintf("  [%d / %d] elapsed=%.1f min, ETA=%.1f min\n",
                w, nrow(windows_tbl), elapsed, eta))
  }
}

# ---- Combine and write ----------------------------------------------------

f3_dt <- rbindlist(results, fill = TRUE)
fwrite(f3_dt, out_csv)

cat("[", format(Sys.time()), "] Wrote:", out_csv,
    "(", nrow(f3_dt), "windows )\n")
