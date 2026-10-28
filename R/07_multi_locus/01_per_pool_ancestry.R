# =============================================================================
# R/07_multi_locus/01_per_pool_ancestry.R
# =============================================================================
# Per-population donor-ancestry test at the three candidate loci. For each
# field pool, alpha is the mean alternate-allele frequency across donor-private
# panel SNPs (pool depth >= 10), and Delta = alpha inside the locus minus alpha
# at panel SNPs outside all three loci. Per locus: one-sided Wilcoxon test of
# median Delta > 0, one-sided regression of Delta on the locus's predictor,
# permutation P values against the width-matched genome-wide null, a joint
# bio05 + pis model with a nested-model F-test, and per-pool 95% bootstrap
# intervals on Delta (B = 1,000).
#
# Manuscript: Figure 5 (drawn by 02_per_pool_ancestry_figure.R); Results
#             "Cx. quinquefasciatus ancestry and sequencing read depth at a CYP6
#             cluster rise with summer temperature", "Cx. quinquefasciatus
#             ancestry is excluded at a developmental-quiescence locus on the
#             same gradient", "A pre-adapted molestus inversion tracks
#             urbanization independently of temperature"; Methods
#             "Per-population donor-ancestry test"; Supplementary Methods
#             "Donor-private ancestry panels and the per-population
#             donor-ancestry test"
#
# Inputs:  results/multilocus/genome_wide_null.rds  (cluster/07_multilocus/01_genome_wide_null.R)
#          data/site_coords_with_spatial_data.csv    (R/00_setup/01_site_predictors.R)
#          deposit: data/vcfs/v4/poolsnp.biallelic.nomiss.vcf.gz
# Outputs: results/multilocus/per_pool_ancestry_results.rds
#          results/multilocus/per_pool_ancestry_per_pool.csv
#            (locus, panel, pool, env, alpha_focal, alpha_genome, delta,
#             delta_ci_lo, delta_ci_hi)
#          results/multilocus/per_pool_ancestry_tests.csv
#            (one row per locus x panel: panel sizes, median Delta, Wilcoxon P,
#             slope and one-sided P, permutation P values, joint-model columns)
#
# Runtime: dominated by reading the whole VCF (11.6 million SNPs x 27 pools)
#          into memory.
# =============================================================================
#
# Notes
# - Focal alpha mirrors the null script: panel SNPs in the region, mean ALT
#   frequency across panel SNPs, depth >= 10. For the non_FL1_deficit panel,
#   alpha = 1 - mean ALT frequency.
# - Bootstrap interval on Delta_i: SNPs are resampled independently in the
#   focal panel and the genome panel, B = 1,000 times; 2.5/97.5% quantiles.
# - Permutation P values use null$<locus>$<primary or secondary> from the null
#   RDS: the fraction of null regions with delta_marginal (or delta_slope) at
#   least as large as the focal value.
# - Focal regions are read from the null RDS and must equal R/helpers/loci.R.

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(VariantAnnotation)
  library(data.table)
})

set.seed(42)


# ---- Configuration ----------------------------------------------------------

VCF_FILE     <- deposit_file('data/vcfs/v4/poolsnp.biallelic.nomiss.vcf.gz')
ENV_FILE     <- 'data/site_coords_with_spatial_data.csv'
NULL_RDS     <- 'results/multilocus/genome_wide_null.rds'
need_files(ENV_FILE, NULL_RDS)

OUT_RDS      <- 'results/multilocus/per_pool_ancestry_results.rds'
OUT_PER_POOL <- 'results/multilocus/per_pool_ancestry_per_pool.csv'
OUT_TESTS    <- 'results/multilocus/per_pool_ancestry_tests.csv'

dir.create(dirname(OUT_RDS), recursive = TRUE, showWarnings = FALSE)

# Bootstrap iterations for the Delta interval
N_BOOT <- 1000L

# Minimum depth for a pool AF to count (matches the null)
MIN_DEPTH <- 10


# =============================================================================
# Stage 1: Load null RDS + extract panel infrastructure
# =============================================================================

cat('[', format(Sys.time()), '] Stage 1: loading null distribution + panel infrastructure\n')

null_data <- readRDS(NULL_RDS)

cfg                <- null_data$config
panel_membership   <- null_data$panel_membership
genome_panel_alpha <- null_data$genome_panel_alpha
field_pools        <- cfg$field_pools
focal_regions      <- cfg$locus_regions

stopifnot(length(field_pools) == 24L)
stopifnot(identical(focal_regions$locus, c('chr1', 'chr2', 'chr3')),
          identical(focal_regions$chrom, c(CHR1, CHR2, CHR3)),
          focal_regions$start == c(CPL1A_START, CPL2A_START, CPL3A_START),
          focal_regions$end   == c(CPL1A_END,   CPL2A_END,   CPL3A_END))
cat(sprintf('  Loaded null with n_regions = %d, %d field pools, %d panels\n',
            cfg$n_random_regions, length(field_pools), length(panel_membership)))
cat('  Focal regions (block scale):\n')
print(focal_regions)


# ---- Locus / donor configuration (mirrors the null script) ------------------

locus_donor_config <- list(
  chr1 = list(donor = 'FL1', panel_name = 'FL1_private', alpha_sign = 1L, env = 'bio05'),
  chr2 = list(donor = 'ME1', panel_name = 'ME1_private', alpha_sign = 1L, env = 'bio05',
              secondary = list(donor = 'FL1', panel_name = 'non_FL1_deficit',
                               alpha_sign = -1L, env = 'bio05',
                               label = 'chr2_nonFL1_deficit')),
  chr3 = list(donor = 'CAL1', panel_name = 'CAL1_private', alpha_sign = 1L, env = 'pis')
)


# =============================================================================
# Stage 2: Read VCF
# =============================================================================

cat('[', format(Sys.time()), '] Stage 2: reading VCF\n')

all_pools <- as.character(samples(scanVcfHeader(VCF_FILE)))

stopifnot(all(REFS %in% all_pools))
stopifnot(setequal(setdiff(all_pools, REFS), field_pools))

field_idx_all <- match(field_pools, all_pools)
ref_idx       <- setNames(match(REFS, all_pools), REFS)

param_all <- ScanVcfParam(samples = all_pools, geno = c('AD','RD','DP'))
vcf_all   <- readVcf(VCF_FILE, param = param_all)

alt_all <- as.matrix(geno(vcf_all)$AD); storage.mode(alt_all) <- 'numeric'
ref_all <- as.matrix(geno(vcf_all)$RD); storage.mode(ref_all) <- 'numeric'
dp_all  <- alt_all + ref_all

gr        <- rowRanges(vcf_all)
snp_chrom <- gsub('|', '_', as.character(seqnames(gr)), fixed = TRUE)
snp_pos   <- as.integer(start(gr))

cat(sprintf('  Loaded %s SNPs x %d pools\n',
            format(nrow(alt_all), big.mark = ','), ncol(alt_all)))

# Panel sizes should match the null
for (p in names(panel_membership)) {
  cat(sprintf('    panel %-18s: %s SNPs (matches: %s)\n',
              p, format(sum(panel_membership[[p]]), big.mark = ','),
              if (sum(panel_membership[[p]]) == null_data$panels[[p]]) 'YES' else 'NO -- VCF may have changed'))
}


# ---- Environmental predictors per field pool --------------------------------

env_dt    <- fread(ENV_FILE)
names(env_dt)[names(env_dt) == 'site'] = 'sample'
env_field <- env_dt[sample %in% field_pools]
env_field <- env_field[match(field_pools, sample)]
stopifnot(identical(env_field$sample, field_pools))


# =============================================================================
# Stage 3: Precompute the focal-exclusion mask (used by genome-panel bootstrap)
# =============================================================================
#
# The same mask the null script used to define in_any_focal, recomputed here
# so it does not have to be stored in the RDS.

in_any_focal <- logical(length(snp_pos))
for (r in seq_len(nrow(focal_regions))) {
  in_any_focal <- in_any_focal |
                  (snp_chrom == focal_regions$chrom[r] &
                   snp_pos >= focal_regions$start[r] &
                   snp_pos <= focal_regions$end[r])
}
cat(sprintf('  SNPs inside any focal region (for genome-panel exclusion): %s\n',
            format(sum(in_any_focal), big.mark = ',')))


# =============================================================================
# Helpers
# =============================================================================

# Compute alpha_focal per pool for a given region + panel.
compute_focal_alpha <- function(chrom, start, end, panel_name, alpha_sign) {
  in_region        <- snp_chrom == chrom & snp_pos >= start & snp_pos <= end
  region_idx       <- which(in_region)
  panel_in_region  <- which(panel_membership[[panel_name]][region_idx])
  panel_region_idx <- region_idx[panel_in_region]

  if (length(panel_region_idx) == 0L) {
    return(list(alpha_per_pool = rep(NA_real_, length(field_idx_all)),
                panel_region_idx = integer(0), n_panel_snps = 0L))
  }

  freq_field <- alt_all[panel_region_idx, field_idx_all, drop = FALSE] /
                dp_all [panel_region_idx, field_idx_all, drop = FALSE]
  freq_field[dp_all[panel_region_idx, field_idx_all, drop = FALSE] < MIN_DEPTH] <- NA_real_
  alpha_per_pool <- colMeans(freq_field, na.rm = TRUE)
  if (alpha_sign == -1L) alpha_per_pool <- 1 - alpha_per_pool

  list(alpha_per_pool   = alpha_per_pool,
       panel_region_idx = panel_region_idx,
       n_panel_snps     = length(panel_region_idx))
}


# Bootstrap CI on Delta_i per pool. Independent resamples in focal and genome
# panels; B replicates per pool.
bootstrap_delta_ci <- function(focal_panel_idx, genome_panel_idx,
                                alpha_sign = 1L, B = N_BOOT, conf = 0.95) {
  n_focal  <- length(focal_panel_idx)
  n_genome <- length(genome_panel_idx)

  if (n_focal == 0L || n_genome == 0L) {
    return(matrix(NA_real_, nrow = length(field_idx_all), ncol = 2,
                  dimnames = list(field_pools, c('lo','hi'))))
  }

  focal_freq <- alt_all[focal_panel_idx, field_idx_all, drop = FALSE] /
                dp_all [focal_panel_idx, field_idx_all, drop = FALSE]
  focal_freq[dp_all[focal_panel_idx, field_idx_all, drop = FALSE] < MIN_DEPTH] <- NA_real_

  genome_freq <- alt_all[genome_panel_idx, field_idx_all, drop = FALSE] /
                 dp_all [genome_panel_idx, field_idx_all, drop = FALSE]
  genome_freq[dp_all[genome_panel_idx, field_idx_all, drop = FALSE] < MIN_DEPTH] <- NA_real_

  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(field_pools))

  for (b in seq_len(B)) {
    f_resample <- sample.int(n_focal,  n_focal,  replace = TRUE)
    g_resample <- sample.int(n_genome, n_genome, replace = TRUE)

    a_focal  <- colMeans(focal_freq [f_resample, , drop = FALSE], na.rm = TRUE)
    a_genome <- colMeans(genome_freq[g_resample, , drop = FALSE], na.rm = TRUE)
    if (alpha_sign == -1L) {
      a_focal  <- 1 - a_focal
      a_genome <- 1 - a_genome
    }

    boot_mat[b, ] <- a_focal - a_genome
  }

  alpha <- (1 - conf) / 2
  ci <- t(apply(boot_mat, 2, function(x) quantile(x, c(alpha, 1 - alpha), na.rm = TRUE)))
  rownames(ci) <- field_pools
  colnames(ci) <- c('lo', 'hi')
  ci
}


# One-sided permutation p (H1: focal > null).
perm_pvalue <- function(focal_stat, null_vec) {
  if (is.na(focal_stat) || length(null_vec) == 0L) return(NA_real_)
  v <- null_vec[is.finite(null_vec)]
  if (length(v) == 0L) return(NA_real_)
  (sum(v >= focal_stat) + 1L) / (length(v) + 1L)
}


# =============================================================================
# Stage 4: Compute focal Delta + tests + permutation p-values per locus
# =============================================================================

cat('[', format(Sys.time()), '] Stage 4: computing focal statistics per locus\n')

per_pool_rows <- list()
test_rows     <- list()

run_one <- function(locus_name, panel_name, alpha_sign, env_var,
                    null_table, locus_label = locus_name) {
  reg <- focal_regions[focal_regions$locus == locus_name, ]

  cat(sprintf('  %s -- panel %-16s -- %s:%d-%d\n',
              locus_label, panel_name, reg$chrom, reg$start, reg$end))

  focal <- compute_focal_alpha(reg$chrom, reg$start, reg$end,
                                panel_name, alpha_sign)
  alpha_genome <- genome_panel_alpha[[panel_name]]

  delta_i <- focal$alpha_per_pool - alpha_genome
  env_vec <- env_field[[env_var]]

  genome_panel_idx <- which(panel_membership[[panel_name]] & !in_any_focal)

  cat(sprintf('    focal panel SNPs: %d   genome panel SNPs: %d\n',
              focal$n_panel_snps, length(genome_panel_idx)))

  ci <- bootstrap_delta_ci(focal$panel_region_idx, genome_panel_idx,
                            alpha_sign = alpha_sign)

  # Per-pool table row
  per_pool_rows[[length(per_pool_rows) + 1L]] <<- data.table(
    locus       = locus_label,
    panel       = panel_name,
    pool        = field_pools,
    env         = env_vec,
    alpha_focal = focal$alpha_per_pool,
    alpha_genome= alpha_genome,
    delta       = delta_i,
    delta_ci_lo = ci[, 'lo'],
    delta_ci_hi = ci[, 'hi']
  )

  # Tests
  delta_clean <- delta_i[is.finite(delta_i)]
  wilcoxon_p <- if (length(delta_clean) >= 3L)
                  wilcox.test(delta_clean, mu = 0, alternative = 'greater')$p.value
                else NA_real_
  median_delta <- median(delta_clean)

  reg_df <- data.frame(delta = delta_i, env = env_vec)
  reg_fit <- tryCatch(lm(delta ~ env, data = reg_df), error = function(e) NULL)
  if (!is.null(reg_fit)) {
    coefs <- summary(reg_fit)$coefficients
    slope_est <- coefs['env', 'Estimate']
    slope_se  <- coefs['env', 'Std. Error']
    slope_p_two <- coefs['env', 'Pr(>|t|)']
    slope_p_one <- if (slope_est > 0) slope_p_two / 2 else 1 - slope_p_two / 2
    slope_ci_lo <- slope_est - 1.96 * slope_se
    slope_ci_hi <- slope_est + 1.96 * slope_se
    adj_r2      <- summary(reg_fit)$adj.r.squared
  } else {
    slope_est <- slope_se <- slope_p_one <- slope_ci_lo <- slope_ci_hi <- adj_r2 <- NA_real_
  }

  perm_p_marginal <- perm_pvalue(median_delta, null_table$delta_marginal)
  perm_p_slope    <- perm_pvalue(slope_est,    null_table$delta_slope)

  # ---- Joint env model -----------------------------------------------------
  # Tests whether the locus's primary predictor effect on Delta survives
  # controlling for the other predictor, and whether adding the other predictor
  # to a primary-only model improves fit (nested-model F-test; columns lrt_*).
  # chr1/chr2 use bio05 as primary and pis as secondary; chr3 uses pis as
  # primary and bio05 as secondary.

  secondary_env_name <- if (env_var == 'bio05') 'pis'
                        else if (env_var == 'pis') 'bio05'
                        else NA_character_
  secondary_env_vec  <- if (!is.na(secondary_env_name)) env_field[[secondary_env_name]] else NULL

  joint_primary_slope      <- NA_real_
  joint_primary_slope_se   <- NA_real_
  joint_primary_slope_p    <- NA_real_
  joint_secondary_slope    <- NA_real_
  joint_secondary_slope_se <- NA_real_
  joint_secondary_slope_p  <- NA_real_
  joint_adj_r2 <- NA_real_
  lrt_F        <- NA_real_
  lrt_p        <- NA_real_

  if (!is.null(secondary_env_vec)) {
    reg_df_joint <- data.frame(delta = delta_i, env = env_vec, env2 = secondary_env_vec)
    # Refit primary-only on the same complete-case set as the joint model so
    # the F-test (anova) is valid (lm's default na.action drops rows with NA
    # in any used variable).
    fit_primary_complete <- tryCatch(lm(delta ~ env,        data = reg_df_joint),
                                      error = function(e) NULL)
    fit_joint            <- tryCatch(lm(delta ~ env + env2, data = reg_df_joint),
                                      error = function(e) NULL)
    if (!is.null(fit_joint)) {
      coefs_joint <- summary(fit_joint)$coefficients
      if ('env'  %in% rownames(coefs_joint)) {
        joint_primary_slope    <- coefs_joint['env',  'Estimate']
        joint_primary_slope_se <- coefs_joint['env',  'Std. Error']
        joint_primary_slope_p  <- coefs_joint['env',  'Pr(>|t|)']
      }
      if ('env2' %in% rownames(coefs_joint)) {
        joint_secondary_slope    <- coefs_joint['env2', 'Estimate']
        joint_secondary_slope_se <- coefs_joint['env2', 'Std. Error']
        joint_secondary_slope_p  <- coefs_joint['env2', 'Pr(>|t|)']
      }
      joint_adj_r2 <- summary(fit_joint)$adj.r.squared
      if (!is.null(fit_primary_complete)) {
        lrt <- tryCatch(anova(fit_primary_complete, fit_joint), error = function(e) NULL)
        if (!is.null(lrt) && nrow(lrt) >= 2L) {
          lrt_F <- lrt$F[2]
          lrt_p <- lrt$`Pr(>F)`[2]
        }
      }
    }
  }

  test_rows[[length(test_rows) + 1L]] <<- data.table(
    locus           = locus_label,
    panel           = panel_name,
    chrom           = reg$chrom,
    start           = reg$start,
    end             = reg$end,
    n_panel_snps    = focal$n_panel_snps,
    n_genome_snps   = length(genome_panel_idx),
    median_delta    = median_delta,
    wilcoxon_p      = wilcoxon_p,
    slope           = slope_est,
    slope_se        = slope_se,
    slope_p_onesided= slope_p_one,
    slope_ci_lo     = slope_ci_lo,
    slope_ci_hi     = slope_ci_hi,
    adj_r2          = adj_r2,
    perm_p_marginal = perm_p_marginal,
    perm_p_slope    = perm_p_slope,
    secondary_env_name       = secondary_env_name,
    joint_primary_slope      = joint_primary_slope,
    joint_primary_slope_se   = joint_primary_slope_se,
    joint_primary_slope_p    = joint_primary_slope_p,
    joint_secondary_slope    = joint_secondary_slope,
    joint_secondary_slope_se = joint_secondary_slope_se,
    joint_secondary_slope_p  = joint_secondary_slope_p,
    joint_adj_r2             = joint_adj_r2,
    lrt_F                    = lrt_F,
    lrt_p                    = lrt_p
  )

  invisible(NULL)
}


for (locus_name in c('chr1', 'chr2', 'chr3')) {
  cfg_l <- locus_donor_config[[locus_name]]

  run_one(locus_name,
          panel_name  = cfg_l$panel_name,
          alpha_sign  = cfg_l$alpha_sign,
          env_var     = cfg_l$env,
          null_table  = null_data$null[[locus_name]]$primary,
          locus_label = locus_name)

  if (!is.null(cfg_l$secondary)) {
    sec <- cfg_l$secondary
    run_one(locus_name,
            panel_name  = sec$panel_name,
            alpha_sign  = sec$alpha_sign,
            env_var     = sec$env,
            null_table  = null_data$null[[locus_name]]$secondary,
            locus_label = sec$label)
  }
}

per_pool_dt <- rbindlist(per_pool_rows)
tests_dt    <- rbindlist(test_rows)


# =============================================================================
# Stage 5: Write tables + save RDS
# =============================================================================

cat('[', format(Sys.time()), '] Stage 5: writing tables\n')

fwrite(per_pool_dt, OUT_PER_POOL)
fwrite(tests_dt,    OUT_TESTS)
cat(sprintf('  Wrote %s (%d rows)\n', OUT_PER_POOL, nrow(per_pool_dt)))
cat(sprintf('  Wrote %s (%d rows)\n', OUT_TESTS,    nrow(tests_dt)))

saveRDS(list(
  config           = list(
    timestamp     = format(Sys.time()),
    n_boot        = N_BOOT,
    min_depth     = MIN_DEPTH,
    null_rds      = NULL_RDS,
    focal_regions = focal_regions
  ),
  per_pool         = per_pool_dt,
  tests            = tests_dt
), OUT_RDS)
cat(sprintf('  Wrote %s\n', OUT_RDS))


# =============================================================================
# Stage 6: Test summary printout
# =============================================================================

cat('\n=== Per-pool ancestry test summary ===\n')
for (i in seq_len(nrow(tests_dt))) {
  r <- tests_dt[i]
  cat(sprintf('\n%s (panel: %s)\n', r$locus, r$panel))
  cat(sprintf('  n_panel_snps : %d  (genome panel: %d)\n', r$n_panel_snps, r$n_genome_snps))
  cat(sprintf('  median Delta : %+.4f\n', r$median_delta))
  cat(sprintf('  Wilcoxon  H1: median > 0 : p = %s\n',
              formatC(r$wilcoxon_p, format = 'g', digits = 3)))
  cat(sprintf('  Regression  Delta ~ env  : slope = %+.5f  (CI [%+.5f, %+.5f]),  adj R^2 = %.3f,  one-sided p = %s\n',
              r$slope, r$slope_ci_lo, r$slope_ci_hi, r$adj_r2,
              formatC(r$slope_p_onesided, format = 'g', digits = 3)))
  cat(sprintf('  Permutation (vs null regions): marginal p = %s,  slope p = %s\n',
              formatC(r$perm_p_marginal, format = 'g', digits = 3),
              formatC(r$perm_p_slope,    format = 'g', digits = 3)))
  if (!is.na(r$secondary_env_name)) {
    primary_env_name <- if (r$secondary_env_name == 'pis') 'bio05'
                        else if (r$secondary_env_name == 'bio05') 'pis'
                        else NA_character_
    cat(sprintf('  Joint  Delta ~ %s + %s  (adj R^2 = %.3f):\n',
                primary_env_name, r$secondary_env_name, r$joint_adj_r2))
    cat(sprintf('    %-5s (partial): slope = %+.5f  (SE %.5f),  p (two-sided) = %s\n',
                primary_env_name,
                r$joint_primary_slope, r$joint_primary_slope_se,
                formatC(r$joint_primary_slope_p, format = 'g', digits = 3)))
    cat(sprintf('    %-5s (partial): slope = %+.5f  (SE %.5f),  p (two-sided) = %s\n',
                r$secondary_env_name,
                r$joint_secondary_slope, r$joint_secondary_slope_se,
                formatC(r$joint_secondary_slope_p, format = 'g', digits = 3)))
    cat(sprintf('    F-test (primary-only vs joint): F = %.2f,  p = %s\n',
                r$lrt_F, formatC(r$lrt_p, format = 'g', digits = 3)))
  }
}

cat('\nDone.\n')
