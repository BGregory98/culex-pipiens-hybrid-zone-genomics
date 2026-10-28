# =============================================================================
# cluster/05_windowed_stats/prepare_pooldata.R
# =============================================================================
# Reads the PoolSNP VCF into a poolfstat pooldata object (27 pools, 80
# haploid genomes each; contig names in underscore form) and builds the table
# of non-overlapping windows of the requested size that hold >= 1,000 SNPs
# (533 windows at 1 Mb).
#
# Step: 7, windowed f3 (cluster/README.md)
# Manuscript: Methods "Windowed f3 admixture scan"
#
# Inputs:  PoolSNP_out_v4/poolsnp.biallelic.nomiss.vcf.gz (step 5)
#          data/metadata/locations_v4.txt (pool-order check only)
#
# Pre-processes the PoolSNP VCF into two RDS files:
#   - pooldata_v4.rds                    (window-size-independent; built once)
#   - windows_tbl_v4_<size>.rds          (built per requested window size)
#
# Splitting these means the expensive VCF parse runs once even when running
# windowed analyses at multiple window sizes (e.g., 1 Mb for the genome scan,
# 250 kb for a focused sensitivity check at the candidate loci).
#
# Usage:
#   module load R
#   # Default: build pooldata (if missing) and windows at 1 Mb
#   Rscript prepare_pooldata.R
#
#   # Build windows at a different size; pooldata is reused if present
#   Rscript prepare_pooldata.R 500000
#   Rscript prepare_pooldata.R 250000
#
#   # Force-rebuild pooldata too (e.g., after re-running PoolSNP):
#   Rscript prepare_pooldata.R 1000000 --rebuild-pooldata
#
# Outputs:
#   SequencingData/data/pooldata_v4.rds
#   SequencingData/data/windows_tbl_v4_<size>.rds
#
# Expected runtime:
#   ~10-15 min for pooldata (one-time, ~11.63M SNPs × 27 pools)
#   ~30 sec for windows_tbl alone (when pooldata is already built)
# =============================================================================

suppressPackageStartupMessages({
  library(poolfstat)
  library(vcfR)
  library(data.table)
  library(dplyr)
  library(tibble)
})

# ---- Parse arguments ------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

WIN_SIZE <- if (length(args) >= 1) as.integer(args[1]) else 1000000L
FORCE_REBUILD_POOLDATA <- "--rebuild-pooldata" %in% args

if (is.na(WIN_SIZE) || WIN_SIZE < 1000) {
  stop("Window size must be a positive integer in bp (e.g., 1000000 for 1 Mb).")
}

# ---- Configuration --------------------------------------------------------

VCF_PATH      <- "PoolSNP_out_v4/poolsnp.biallelic.nomiss.vcf.gz"
SAMPLES_PATH  <- "data/metadata/locations_v4.txt"
OUT_POOLDATA  <- "SequencingData/data/pooldata_v4.rds"

MIN_SNPS_PER_WIN <- 1000L
HAPLOID_SIZE     <- 80L

# Format window size for filename: 1000000 -> "1M", 500000 -> "500K", etc.
format_win_size <- function(n) {
  if (n %% 1e6 == 0)  return(paste0(n %/% 1e6, "M"))
  if (n %% 1e3 == 0)  return(paste0(n %/% 1e3, "K"))
  return(as.character(n))
}
WIN_SIZE_LABEL <- format_win_size(WIN_SIZE)
OUT_WINDOWS    <- sprintf("SequencingData/data/windows_tbl_v4_%s.rds", WIN_SIZE_LABEL)

cat("Configuration:\n")
cat("  Window size:  ", WIN_SIZE, "bp (", WIN_SIZE_LABEL, ")\n")
cat("  Min SNPs/win: ", MIN_SNPS_PER_WIN, "\n")
cat("  Pooldata RDS: ", OUT_POOLDATA, "\n")
cat("  Windows RDS:  ", OUT_WINDOWS, "\n\n")

# ---- Canonicalize chromosome names to underscore form -------------------

canonicalize_chr <- function(x) {
  ifelse(grepl("^ENA\\|", x), gsub("|", "_", x, fixed = TRUE), x)
}

# ==== Stage 1: Build pooldata RDS (if missing or forced) =================

if (file.exists(OUT_POOLDATA) && !FORCE_REBUILD_POOLDATA) {
  cat("[", format(Sys.time()), "] Loading existing pooldata from",
      OUT_POOLDATA, "...\n")
  pooldata <- readRDS(OUT_POOLDATA)
  cat("  pooldata: ", pooldata@nsnp, "SNPs x", pooldata@npools, "pools\n\n")
} else {

  if (!file.exists(VCF_PATH)) {
    stop("VCF not found at ", VCF_PATH,
         ". Has v4 PoolSNP finished? Are you in the project root?")
  }

  cat("[", format(Sys.time()), "] Reading VCF header to get pool names...\n")
  vcf_header <- read.vcfR(VCF_PATH, verbose = FALSE)
  pool_names <- colnames(vcf_header@gt)[-1]
  n_pools    <- length(pool_names)
  cat("  Pools in VCF (n =", n_pools, "):\n   ",
      paste(pool_names, collapse = ", "), "\n")

  if (file.exists(SAMPLES_PATH)) {
    expected <- readLines(SAMPLES_PATH) |> trimws()
    expected <- expected[nzchar(expected)]
    if (!identical(pool_names, expected)) {
      cat("  WARNING: VCF pool order does not match samples.txt order.\n")
      cat("    Using VCF order as authoritative for pooldata indexing.\n")
    } else {
      cat("  Pool order matches samples.txt. ✓\n")
    }
  }
  rm(vcf_header); gc()

  cat("[", format(Sys.time()),
      "] Parsing VCF into pooldata (this is the slow step)...\n")
  pooldata <- vcf2pooldata(
    vcf.file  = VCF_PATH,
    poolsizes = rep(HAPLOID_SIZE, n_pools),
    poolnames = pool_names,
    verbose   = FALSE
  )
  cat("  pooldata: ", pooldata@nsnp, "SNPs x", pooldata@npools, "pools\n")

  orig_chroms <- unique(pooldata@snp.info$Chromosome)
  cat("  Chromosome names in VCF:", paste(orig_chroms, collapse = ", "), "\n")
  pooldata@snp.info$Chromosome <- canonicalize_chr(pooldata@snp.info$Chromosome)
  new_chroms <- unique(pooldata@snp.info$Chromosome)
  if (!identical(orig_chroms, new_chroms)) {
    cat("  Canonicalized to:", paste(new_chroms, collapse = ", "), "\n")
  }

  saveRDS(pooldata, OUT_POOLDATA)
  cat("  Wrote:", OUT_POOLDATA,
      "(", round(file.size(OUT_POOLDATA) / 1e6, 1), "MB )\n\n")
}

# ==== Stage 2: Build windows_tbl for the requested size ==================

cat("[", format(Sys.time()), "] Building", WIN_SIZE_LABEL,
    "windows (≥", MIN_SNPS_PER_WIN, "SNPs)...\n")

snp_info <- as.data.frame(pooldata@snp.info)
colnames(snp_info) <- tolower(colnames(snp_info))
stopifnot(all(c("chromosome", "position") %in% colnames(snp_info)))

snp_tbl <- tibble(
  snp_index  = seq_len(nrow(snp_info)),
  chromosome = as.factor(snp_info$chromosome),
  position   = as.integer(snp_info$position)
)

windows_tbl <- snp_tbl %>%
  mutate(bin = (position - 1L) %/% WIN_SIZE) %>%
  group_by(chromosome, bin) %>%
  reframe(
    snp_idx = list(snp_index),
    n_snps  = dplyr::n(),
    start   = bin[1] * WIN_SIZE + 1L,
    end     = (bin[1] + 1L) * WIN_SIZE
  ) %>%
  dplyr::filter(n_snps >= MIN_SNPS_PER_WIN) %>%
  arrange(chromosome, start)

cat("  Windows passing min-SNPs filter:", nrow(windows_tbl), "\n")
cat("  Per chromosome:\n")
windows_tbl %>%
  dplyr::count(chromosome) %>%
  rowwise() %>%
  group_walk(~ cat("    ", as.character(.x$chromosome), ":", .x$n, "windows\n"))

saveRDS(windows_tbl, OUT_WINDOWS)
cat("  Wrote:", OUT_WINDOWS,
    "(", round(file.size(OUT_WINDOWS) / 1e3, 1), "KB )\n\n")

cat("[", format(Sys.time()), "] Done.\n")
