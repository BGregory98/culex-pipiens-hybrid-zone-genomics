# =============================================================================
# R/07_multi_locus/03_haplotype_block_diagnostics.R
# =============================================================================
# Haplotype-block diagnostics at the three candidate loci, each tested against
# the genome-wide null of width-matched random regions:
#   (1) pc1_pct_var   : PC1 % variance of the field-pool AF matrix at region SNPs
#   (2) pc1_fst_cor   : |cor(PC1 score, mean FST(pool, donor))| across pools
#   (3) fst_donor_var : variance of mean FST(pool, donor) across the 24 pools
# with one-sided upper-tail permutation P values read from the null RDS.
#
# Manuscript: Results "A pre-adapted molestus inversion tracks urbanization
#             independently of temperature" (PC1 %, |r|); Methods "Locus-level
#             diagnostics against a genome-wide null"; Supplementary Methods
#             "Haplotype-block diagnostics"; Supplementary Results
#             "Haplotype-block diagnostics at the three candidate loci"
#
# Inputs:  results/multilocus/genome_wide_null.rds  (cluster/07_multilocus/01_genome_wide_null.R)
#          deposit: data/vcfs/v4/poolsnp.biallelic.nomiss.vcf.gz (+ .tbi)
# Outputs: results/multilocus/multi_locus_pca_fst_summary.csv
#          results/multilocus/multi_locus_pca_fst_scores.csv
#
# Runtime: minutes (region-restricted VCF reads)
# =============================================================================
#
# Three tested rows (chr2 is read on the QUI-versus-non-QUI axis):
#   chr1 / FL1            (null chr1 primary)
#   chr2 / FL1 (deficit)  (null chr2 secondary)
#   chr3 / CAL1           (null chr3 primary)
#
# The FST estimator and the PCA filters must match the null script
# (cluster/07_multilocus/01_genome_wide_null.R) exactly for the permutation
# P values to be valid: per-SNP Reynolds-style FST averaged across SNPs, and
# PCA on SNPs with field-pool mean frequency in [0.05, 0.95] and no masked
# (depth < 10) cells.

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(VariantAnnotation)
})


# ---- Configuration ----------------------------------------------------------

VCF_FILE  <- deposit_file('data/vcfs/v4/poolsnp.biallelic.nomiss.vcf.gz')
NULL_RDS  <- 'results/multilocus/genome_wide_null.rds'
need_files(NULL_RDS)

SUMMARY_OUT <- 'results/multilocus/multi_locus_pca_fst_summary.csv'
SCORES_OUT  <- 'results/multilocus/multi_locus_pca_fst_scores.csv'

dir.create(dirname(SUMMARY_OUT), recursive = TRUE, showWarnings = FALSE)

# Filters matching the null
MAF_MIN   <- 0.05
MIN_DEPTH <- 10


# ---- Locus / panel test configuration ---------------------------------------

TEST_CFG <- list(
  list(label    = 'chr1',
       locus    = 'chr1',
       chrom    = CHR1,
       start    = CPL1A_START, end = CPL1A_END,
       donor    = 'FL1',
       null_key = c('chr1','primary')),
  list(label    = 'chr2_nonFL1_deficit',
       locus    = 'chr2',
       chrom    = CHR2,
       start    = CPL2A_START, end = CPL2A_END,
       donor    = 'FL1',
       null_key = c('chr2','secondary')),
  list(label    = 'chr3',
       locus    = 'chr3',
       chrom    = CHR3,
       start    = CPL3A_START, end = CPL3A_END,
       donor    = 'CAL1',
       null_key = c('chr3','primary'))
)


# =============================================================================
# Helpers (identical to the null script where they overlap)
# =============================================================================

# Per-SNP Reynolds-style FST to donor pool, averaged across SNPs per field pool.
# Same estimator as the null's compute_mean_fst_to_donor; do not modify
# without also updating the null and rerunning it.
compute_mean_fst_to_donor <- function(alt, dp, donor_idx, field_idx,
                                      min_depth = MIN_DEPTH) {
  freq <- alt / dp
  freq[dp < min_depth] <- NA_real_
  p_donor <- freq[, donor_idx]
  vapply(field_idx, function(j) {
    p_pool <- freq[, j]
    pbar   <- (p_pool + p_donor) / 2
    num    <- (p_pool - p_donor)^2
    den    <- pbar * (1 - pbar)
    valid  <- is.finite(num) & is.finite(den) & den > 0
    if (sum(valid) == 0L) return(NA_real_)
    mean(num[valid] / den[valid])
  }, numeric(1))
}

# PCA on per-pool AF matrix at field pools, MAF + complete-rows filtered.
# Same filters as the null's compute_region_pca; also returns PC2.
compute_region_pca <- function(alt, dp, field_idx, min_depth = MIN_DEPTH,
                               maf_min = MAF_MIN) {
  freq <- alt[, field_idx] / dp[, field_idx]
  freq[dp[, field_idx] < min_depth] <- NA_real_

  mean_freq <- rowMeans(freq, na.rm = TRUE)
  keep_maf  <- !is.na(mean_freq) & mean_freq >= maf_min & mean_freq <= 1 - maf_min
  freq      <- freq[keep_maf, , drop = FALSE]
  complete  <- rowSums(is.na(freq)) == 0
  freq      <- freq[complete, , drop = FALSE]

  if (nrow(freq) < 10L || ncol(freq) < 3L) {
    return(list(pct_var_PC1 = NA_real_,
                pct_var_PC2 = NA_real_,
                scores_PC1  = rep(NA_real_, length(field_idx)),
                scores_PC2  = rep(NA_real_, length(field_idx)),
                n_snps      = nrow(freq)))
  }

  centered <- freq - rowMeans(freq)
  pca      <- prcomp(t(centered), center = FALSE, scale. = FALSE)
  pct_var  <- 100 * (pca$sdev^2) / sum(pca$sdev^2)

  list(pct_var_PC1 = pct_var[1],
       pct_var_PC2 = pct_var[2],
       scores_PC1  = pca$x[, 1],
       scores_PC2  = pca$x[, 2],
       n_snps      = nrow(freq))
}

# Read AD/RD matrices from a region of the v4 VCF for the listed pools.
# Returns alt and dp matrices (rows = SNPs in region, cols = pools). The VCF
# contigs use the underscore names, so the region is queried as given.
read_region_AD_RD <- function(vcf_path, chrom, start, end, pool_names) {
  region_gr <- GRanges(seqnames = chrom,
                       ranges   = IRanges(start = start, end = end))
  param <- ScanVcfParam(samples = pool_names,
                        geno    = c('AD', 'RD'),
                        which   = region_gr)
  vcf <- readVcf(vcf_path, param = param)
  alt <- as.matrix(geno(vcf)$AD); storage.mode(alt) <- 'numeric'
  ref <- as.matrix(geno(vcf)$RD); storage.mode(ref) <- 'numeric'
  dp  <- alt + ref
  list(alt = alt, dp = dp, n_snps_raw = nrow(alt))
}

# Upper-tailed permutation p-value: probability that null statistic >= focal.
# Uses the (1 + #ge) / (1 + n) form to avoid p = 0.
perm_p_upper <- function(focal, null_vec) {
  null_vec <- null_vec[is.finite(null_vec)]
  if (length(null_vec) == 0L || !is.finite(focal)) return(NA_real_)
  (1L + sum(null_vec >= focal)) / (1L + length(null_vec))
}


# =============================================================================
# Stage 1: Load the null + extract config
# =============================================================================

cat('[', format(Sys.time()), '] Stage 1: loading genome-wide null distribution\n')
null_obj <- readRDS(NULL_RDS)

field_pools <- null_obj$config$field_pools
stopifnot(length(field_pools) == 24L)
cat(sprintf('  Field pools (n=%d): %s ...\n', length(field_pools),
            paste(head(field_pools, 5), collapse = ', ')))

# Field pools plus the references (the donor is one of the references)
ALL_POOLS <- c(field_pools, REFS)


# =============================================================================
# Stage 2: For each test row, read region, compute statistics, look up perm p
# =============================================================================

cat('[', format(Sys.time()), '] Stage 2: computing focal statistics per locus\n')

summary_rows <- list()
scores_rows  <- list()

for (cfg in TEST_CFG) {

  cat(sprintf('\n--- %s (donor = %s; null = %s/%s) ---\n',
              cfg$label, cfg$donor, cfg$null_key[1], cfg$null_key[2]))
  cat(sprintf('    %s:%.2f-%.2f Mb\n', cfg$chrom,
              cfg$start/1e6, cfg$end/1e6))

  # Read region from VCF (region-restricted).
  reg <- read_region_AD_RD(VCF_FILE, cfg$chrom, cfg$start, cfg$end,
                           pool_names = ALL_POOLS)
  cat(sprintf('    %s SNPs in region (raw)\n',
              format(reg$n_snps_raw, big.mark = ',')))

  # Column indices into the region matrices (readVcf returns columns in the
  # order requested via `samples`).
  field_idx_local <- match(field_pools, ALL_POOLS)
  donor_idx_local <- match(cfg$donor,    ALL_POOLS)
  stopifnot(!any(is.na(field_idx_local)), !is.na(donor_idx_local))
  stopifnot(identical(colnames(reg$alt), ALL_POOLS))

  pca <- compute_region_pca(reg$alt, reg$dp, field_idx = field_idx_local)

  fst_per_pool <- compute_mean_fst_to_donor(reg$alt, reg$dp,
                                            donor_idx = donor_idx_local,
                                            field_idx = field_idx_local)

  pc1_fst_cor_signed <- if (all(is.finite(pca$scores_PC1)) &&
                            all(is.finite(fst_per_pool)))
                          cor(pca$scores_PC1, fst_per_pool) else NA_real_
  pc1_fst_cor_abs    <- if (is.na(pc1_fst_cor_signed)) NA_real_
                        else abs(pc1_fst_cor_signed)
  fst_donor_var      <- var(fst_per_pool, na.rm = TRUE)

  # Permutation p-values from the null
  null_dt <- null_obj$null[[ cfg$null_key[1] ]][[ cfg$null_key[2] ]]

  perm_p_pc1_var <- perm_p_upper(pca$pct_var_PC1,  null_dt$pc1_pct_var)
  perm_p_pc1_fst <- perm_p_upper(pc1_fst_cor_abs,  null_dt$pc1_fst_cor)
  perm_p_fst_var <- perm_p_upper(fst_donor_var,    null_dt$fst_donor_var)

  n_null_valid <- sum(is.finite(null_dt$pc1_pct_var))

  cat(sprintf('    n_snps (post-filter)    : %d\n', pca$n_snps))
  cat(sprintf('    n_null_regions_valid    : %d\n', n_null_valid))
  cat(sprintf('    PC1 %% var               : %.2f  (perm p = %.4g)\n',
              pca$pct_var_PC1, perm_p_pc1_var))
  cat(sprintf('    |PC1-FST(%-4s) cor|      : %.3f (perm p = %.4g)\n',
              cfg$donor, pc1_fst_cor_abs, perm_p_pc1_fst))
  cat(sprintf('    FST(pool, %-4s) variance : %.5f (perm p = %.4g)\n',
              cfg$donor, fst_donor_var, perm_p_fst_var))

  summary_rows[[cfg$label]] <- tibble(
    label              = cfg$label,
    locus              = cfg$locus,
    donor              = cfg$donor,
    chrom              = cfg$chrom,
    start              = cfg$start,
    end                = cfg$end,
    n_snps             = pca$n_snps,
    pc1_pct_var        = pca$pct_var_PC1,
    pc2_pct_var        = pca$pct_var_PC2,
    pc1_fst_cor_signed = pc1_fst_cor_signed,
    pc1_fst_cor_abs    = pc1_fst_cor_abs,
    fst_donor_var      = fst_donor_var,
    fst_donor_mean     = mean(fst_per_pool, na.rm = TRUE),
    n_null_regions     = n_null_valid,
    perm_p_pc1_var     = perm_p_pc1_var,
    perm_p_pc1_fst     = perm_p_pc1_fst,
    perm_p_fst_var     = perm_p_fst_var
  )

  scores_rows[[cfg$label]] <- tibble(
    label    = cfg$label,
    locus    = cfg$locus,
    donor    = cfg$donor,
    pool     = field_pools,
    PC1      = pca$scores_PC1,
    PC2      = pca$scores_PC2,
    mean_fst = fst_per_pool / 4   # 4*GST -> Nei's pairwise GST
  )
}


# =============================================================================
# Stage 3: Combined outputs
# =============================================================================

summary_tbl <- bind_rows(summary_rows)
scores_tbl  <- bind_rows(scores_rows)

cat('\n=== COMBINED SUMMARY ===\n')
print(as.data.frame(summary_tbl %>%
                      dplyr::select(label, n_snps, pc1_pct_var, pc1_fst_cor_abs,
                             fst_donor_var, perm_p_pc1_var, perm_p_pc1_fst,
                             perm_p_fst_var)),
      row.names = FALSE, digits = 4)

write.csv(summary_tbl, SUMMARY_OUT, row.names = FALSE)
write.csv(scores_tbl,  SCORES_OUT,  row.names = FALSE)
cat(sprintf('\nWrote %s\n', SUMMARY_OUT))
cat(sprintf('Wrote %s\n',   SCORES_OUT))

cat(sprintf('\n[%s] Done.\n', format(Sys.time())))
