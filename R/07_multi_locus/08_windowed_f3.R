# =============================================================================
# R/07_multi_locus/08_windowed_f3.R
# =============================================================================
# Analyses of the windowed f3 statistics (1 Mb windows; each field pool tested
# against the three reference pairings; BH-adjusted admixture calls made on the
# cluster):
#   (A) Genome-wide: per-population admixed-window counts per reference pair,
#       regressed on each environmental predictor, and joint bio05 + pis
#       models with a nested-model F-test.
#   (B) Per focal cluster (CpL1a = cluster1, CpL2a = cluster2): each focal
#       1 Mb window's z-score against the same pool's genome-wide background
#       (empirical percentile), and the focal-minus-background-median deviation
#       regressed on the cluster's predictor (bio05) per reference pair.
# Site predictors are the 500 m buffer means of Table S1. Reference pairs:
# FL1-CAL1 = QUI-MOL, ME1-CAL1 = PIP-MOL, ME1-FL1 = PIP-QUI.
#
# Manuscript: Results "Genome-wide admixture is structured by both temperature
#             and urbanization" (533 windows, admixed-window means and
#             regressions); Results "Cx. quinquefasciatus ancestry and
#             sequencing read depth at a CYP6 cluster rise with summer
#             temperature" (CpL1a QUI-MOL beta); Methods "Windowed f3 admixture
#             scan"; Supplementary Results "CpL1a is uniquely enriched for
#             QUI--MOL admixture" and "The f3 topology at CpL2a differs from
#             CpL1a"; Figure S10 (drawn by 09_windowed_f3_figure.R)
#
# Inputs:  data/windowed_F3_all_comparisons_v4_1M.csv  (cluster/05_windowed_stats/combine_windowed_outputs.R)
#          data/site_coords_with_spatial_data.csv       (R/00_setup/01_site_predictors.R)
#          data/sample_categories.csv                   (via R/helpers/sample_metadata.R)
# Outputs: results/f3/focal_window_f3_cluster1.csv, results/f3/focal_window_f3_cluster2.csv
#            (central focal window per pool and reference pair; read by
#             09_windowed_f3_figure.R)
#
# Runtime: about a minute
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/cluster_presets.R")
source("R/helpers/sample_metadata.R")

# ---- Libraries --------------------------------------------------------------

library(tidyverse)
library(sf)

# Focal clusters the paper reports (CpL1a, CpL2a)
CLUSTERS <- c('cluster1', 'cluster2')

F3_FILE  <- 'data/windowed_F3_all_comparisons_v4_1M.csv'
ENV_FILE <- 'data/site_coords_with_spatial_data.csv'
need_files(F3_FILE, ENV_FILE)

dir.create('results/f3', recursive = TRUE, showWarnings = FALSE)

# The 1 Mb f3 windows that overlap a cluster. Windows are indexed by their
# 1-based start position (X000001), matching the windowed f3 output
# convention. Any window that overlaps [cluster_start, cluster_end] is retained.
get_cluster_window_starts <- function(start, end, window_size = 1e6) {
  first <- floor(start / window_size) * window_size + 1
  last  <- floor((end - 1) / window_size) * window_size + 1
  seq(first, last, by = window_size)
}


# ---- Sample and environmental data ------------------------------------------

# Site-level predictors: means within a 500 m buffer around each site
# (Table S1).
site_coords <- read.csv(ENV_FILE) %>%
  dplyr::filter(site %in% field_names) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)

env_data <- site_coords %>%
  as.data.frame() %>%
  dplyr::rename(popA = site)


# ---- Load windowed f3 results -----------------------------------------------

# Windowed f3 statistics per population for the three reference pairs.
f3_data_frame <- read.csv(F3_FILE)
f3_data_frame$ref = as.factor(paste0(f3_data_frame$popB, '-', f3_data_frame$popC))

# Genome-wide factor ordering of pools by bio05 (coolest first).
bio05_order <- env_data %>%
  arrange(bio05) %>%
  pull(popA)

f3_data_frame$popA = factor(f3_data_frame$popA, levels = bio05_order)

# How many 1 Mb windows were retained genome-wide?
n_windows_retained <- f3_data_frame %>%
  distinct(window_chr, window_start) %>%
  nrow()


# =============================================================================
# (A) Per-population admixed-window counts vs. environmental variables
# =============================================================================
#
# If MOL introgression is structured by urbanization rather than (or in
# addition to) temperature, the number of admixed windows per population in
# MOL-containing contrasts (ME1-CAL1, FL1-CAL1) should correlate with
# urbanization metrics, while ME1-FL1 (the MOL-free contrast) should not.
# =============================================================================

# ---- Count admixed windows per population per reference pair ----------------

admixed_counts <- f3_data_frame %>%
  dplyr::filter(popA %in% field_names) %>%
  dplyr::rename(sample = popA) %>%
  group_by(ref, sample) %>%
  summarize(n_admixed_windows = sum(admix_bh, na.rm = TRUE), .groups = 'drop') %>%
  left_join(sample_info_all, by = 'sample')

# Summary statistics per reference pair.
cat('\n=== Admixed-window counts per population, by reference pair ===\n')
admixed_summary <- admixed_counts %>%
  group_by(ref) %>%
  summarize(
    mean_n   = mean(n_admixed_windows),
    sd_n     = sd(n_admixed_windows),
    min_n    = min(n_admixed_windows),
    max_n    = max(n_admixed_windows),
    total_n  = sum(n_admixed_windows),
    .groups  = 'drop'
  )
print(admixed_summary)

# ---- Robustness to the grouping of the BH correction -------------------------
# The stored calls correct within each reference pairing (all 24 populations as
# one family). Recompute the per-population mean admixed-window count under a
# single correction across all tests, and under one correction per population x
# pairing (Supplementary Results, robustness of the windowed f3 correction).
field_tests <- f3_data_frame %>% dplyr::filter(popA %in% field_names)
bh_mean <- function(d) {
  d %>%
    mutate(admix = p_adj < 0.05) %>%
    group_by(ref, popA) %>% summarize(n = sum(admix, na.rm = TRUE), .groups = 'drop') %>%
    group_by(ref) %>% summarize(mean_n = mean(n), .groups = 'drop')
}
bh_robust <- bind_rows(
  field_tests %>% mutate(p_adj = p.adjust(p_one_sided, method = 'BH')) %>%
    bh_mean() %>% mutate(scheme = 'single correction, all tests'),
  field_tests %>% group_by(ref) %>% mutate(p_adj = p.adjust(p_one_sided, method = 'BH')) %>%
    ungroup() %>% bh_mean() %>% mutate(scheme = 'per reference pairing (as used)'),
  field_tests %>% group_by(ref, popA) %>% mutate(p_adj = p.adjust(p_one_sided, method = 'BH')) %>%
    ungroup() %>% bh_mean() %>% mutate(scheme = 'per population x pairing')
) %>%
  tidyr::pivot_wider(names_from = ref, values_from = mean_n)
cat('\n=== Mean admixed windows per population under three BH groupings ===\n')
print(as.data.frame(bh_robust), digits = 4)


# ---- Join to environmental data ---------------------------------------------

env_flat <- env_data %>%
  dplyr::rename(sample = popA) %>%
  dplyr::select(-geometry) %>%
  as.data.frame()

counts_env <- admixed_counts %>%
  left_join(env_flat, by = 'sample')


# ---- Univariate environmental regressions per reference pair ----------------
#
# For each reference pair, regress per-population admixed-window count on each
# environmental predictor separately and report slope, adjusted R-squared,
# and p-value.

env_predictors <- c('bio05', 'bio07', 'bio18',
                    'pis', 'tree_canopy', 'pd',
                    'elevation', 'wtd')

fit_univariate_env <- function(dat, predictors) {
  lapply(predictors, function(pred) {
    if (!pred %in% names(dat)) return(NULL)
    fit <- lm(reformulate(pred, response = 'n_admixed_windows'), data = dat)
    sm  <- summary(fit)
    co  <- coef(sm)
    slope <- if (pred %in% rownames(co)) co[pred, 'Estimate']   else NA_real_
    p_val <- if (pred %in% rownames(co)) co[pred, 'Pr(>|t|)']   else NA_real_
    data.frame(
      predictor = pred,
      slope     = slope,
      p_val     = p_val,
      adj_r2    = sm$adj.r.squared,
      n         = nrow(dat),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
}

cat('\n=== Univariate environmental predictors of admixed-window count ===\n')
univariate_results <- lapply(levels(counts_env$ref), function(r) {
  dat <- counts_env %>% dplyr::filter(ref == r)
  fit_univariate_env(dat, env_predictors) %>%
    mutate(ref = r) %>%
    arrange(desc(adj_r2))
}) %>% setNames(levels(counts_env$ref))

for (r in names(univariate_results)) {
  cat('\n---', r, '---\n')
  print(univariate_results[[r]])
}


# ---- Joint (multivariate) models: is urbanization effect independent? ------
#
# For each reference pair, fit two nested models:
#   Model 1: n_admixed_windows ~ bio05
#   Model 2: n_admixed_windows ~ bio05 + pis
# A nested-model F-test assesses whether adding impervious surface improves
# fit over temperature alone.

cat('\n=== Joint models: bio05 vs bio05 + pis ===\n')
joint_results <- list()
for (r in c('FL1-CAL1', 'ME1-CAL1', 'ME1-FL1')) {
  dat <- counts_env %>% dplyr::filter(ref == r)

  m_bio05      <- lm(n_admixed_windows ~ bio05,       data = dat)
  m_bio05_pis  <- lm(n_admixed_windows ~ bio05 + pis, data = dat)
  m_pis        <- lm(n_admixed_windows ~ pis,         data = dat)

  cat('\n---', r, '---\n')
  cat('Model: bio05 only\n');    print(summary(m_bio05)$coefficients);     cat('Adj R^2 =', summary(m_bio05)$adj.r.squared, '\n')
  cat('\nModel: pis only\n');    print(summary(m_pis)$coefficients);       cat('Adj R^2 =', summary(m_pis)$adj.r.squared, '\n')
  cat('\nModel: bio05 + pis\n'); print(summary(m_bio05_pis)$coefficients); cat('Adj R^2 =', summary(m_bio05_pis)$adj.r.squared, '\n')
  cat('\nF-test: does adding pis improve over bio05 alone?\n')
  print(anova(m_bio05, m_bio05_pis))

  joint_results[[r]] <- list(bio05 = m_bio05, pis = m_pis, both = m_bio05_pis)
}


# =============================================================================
# MANUSCRIPT VALUES (A): genome-wide admixture
# =============================================================================

cat('\n\n')
cat('============================================================\n')
cat('MANUSCRIPT VALUES (A): genome-wide admixture\n')
cat('============================================================\n')

cat(sprintf('\nRetained 1 Mb windows genome-wide: %d\n', n_windows_retained))
cat(sprintf('Total tests (windows x pops x 3 ref pairs): %d\n',
            nrow(f3_data_frame %>% dplyr::filter(popA %in% field_names))))

cat('\n--- Admixed-window counts per reference pair (BH < 0.05) ---\n')
print(admixed_summary)

cat('\n--- Univariate predictors of admixed-window count, per reference pair ---\n')
for (r in names(univariate_results)) {
  cat('\n  >>>', r, '<<<\n')
  print(univariate_results[[r]], row.names = FALSE)
}

cat('\n--- Joint models (bio05 vs bio05 + pis), per reference pair ---\n')
for (r in c('FL1-CAL1', 'ME1-CAL1', 'ME1-FL1')) {
  fits <- joint_results[[r]]
  sm_bio05     <- summary(fits$bio05)
  sm_pis       <- summary(fits$pis)
  sm_both      <- summary(fits$both)
  lrt          <- anova(fits$bio05, fits$both)
  cat(sprintf('\n  >>> %s <<<\n', r))
  cat(sprintf('    bio05 alone:  beta = %.3f, adj R^2 = %.3f, p = %.3g\n',
              coef(sm_bio05)['bio05','Estimate'],
              sm_bio05$adj.r.squared,
              coef(sm_bio05)['bio05','Pr(>|t|)']))
  cat(sprintf('    pis   alone:  beta = %.3f, adj R^2 = %.3f, p = %.3g\n',
              coef(sm_pis)['pis','Estimate'],
              sm_pis$adj.r.squared,
              coef(sm_pis)['pis','Pr(>|t|)']))
  cat(sprintf('    joint model:  bio05 beta = %.3f (p = %.3g), pis beta = %.3f (p = %.3g), adj R^2 = %.3f\n',
              coef(sm_both)['bio05','Estimate'],
              coef(sm_both)['bio05','Pr(>|t|)'],
              coef(sm_both)['pis','Estimate'],
              coef(sm_both)['pis','Pr(>|t|)'],
              sm_both$adj.r.squared))
  cat(sprintf('    F-test (add pis over bio05 alone): F = %.3f, p = %.3g\n',
              lrt$F[2], lrt$`Pr(>F)`[2]))
}

cat('\n============================================================\n')
cat('END MANUSCRIPT VALUES (A)\n')
cat('============================================================\n\n')


# =============================================================================
# (B) Focal-cluster candidate window analysis, one cluster at a time
# =============================================================================

# For a given focal window, compare its f3 z-score to the empirical distribution
# of z-scores across all other windows (same population x reference pair),
# returning focal value, background summary statistics, and empirical p-value.
run_empirical <- function(f3, chrom_focal, focal_start) {
  f3 %>%
    group_by(ref, popA) %>%
    group_modify(~{
      focal <- .x %>% dplyr::filter(window_chr == chrom_focal, window_start == focal_start)
      bg    <- .x %>% dplyr::filter(!(window_chr == chrom_focal & window_start == focal_start))

      tibble(
        focal_start    = focal_start,
        focal_z        = focal$z_score[1],
        focal_est      = focal$estimate[1],
        focal_admix    = focal$admix_bh[1],
        bg_mean_z      = mean(bg$z_score, na.rm = TRUE),
        bg_median_z    = median(bg$z_score, na.rm = TRUE),
        bg_q025_z      = quantile(bg$z_score, 0.025, na.rm = TRUE),
        bg_q975_z      = quantile(bg$z_score, 0.975, na.rm = TRUE),
        bg_admix_rate  = mean(bg$admix_bh, na.rm = TRUE),
        # Empirical one-sided p-value: proportion of background windows
        # with z-score at least as low as focal (with +1 continuity correction).
        emp_p_z        = (1 + sum(bg$z_score <= focal$z_score[1], na.rm = TRUE)) / (nrow(bg) + 1),
        emp_percentile = 1 - (1 + sum(bg$z_score <= focal$z_score[1], na.rm = TRUE)) / (nrow(bg) + 1),
        emp_rank       = 1 + sum(bg$z_score <= focal$z_score[1], na.rm = TRUE),
        n_bg           = nrow(bg)
      )
    }) %>%
    ungroup()
}

run_cluster <- function(active_cluster) {

  cfg                <- cluster_presets[[active_cluster]]
  cluster_chr        <- cfg$chr
  cluster_start      <- cfg$start
  cluster_end        <- cfg$end
  cluster_predictor  <- cfg$predictor
  cluster_pred_label <- cfg$predictor_label
  file_suffix        <- paste0('_', active_cluster)

  cluster_windows <- get_cluster_window_starts(cluster_start, cluster_end)

  # The central focal window is the middle overlapping window of the cluster.
  focal_window_for_plot <- cluster_windows[ceiling(length(cluster_windows) / 2)]

  cat(sprintf('\nActive cluster: %s  (%s, %.2f-%.2f Mb)\n',
              active_cluster, cluster_chr,
              cluster_start / 1e6, cluster_end / 1e6))
  cat(sprintf('Predictor:      %s  (%s)\n', cluster_predictor, cluster_pred_label))
  cat(sprintf('Overlapping 1 Mb f3 windows (n = %d): %s\n',
              length(cluster_windows),
              paste(cluster_windows, collapse = ', ')))
  cat(sprintf('Central focal window: %d\n\n', focal_window_for_plot))

  # Cluster-predictor order of pools (ascending).
  cluster_pred_order <- env_data %>%
    arrange(.data[[cluster_predictor]]) %>%
    pull(popA)

  # Restrict to field populations, ordered by the cluster's predictor.
  f3 <- f3_data_frame %>%
    dplyr::filter(popA %in% field_names) %>%
    mutate(popA = factor(as.character(popA), levels = cluster_pred_order))

  # Run `run_empirical` once per overlapping 1 Mb window in the focal cluster,
  # then attach environmental data for the regressions.
  empirical_results <- bind_rows(
    lapply(cluster_windows, function(ws) run_empirical(f3, cluster_chr, ws))
  ) %>%
    left_join(env_data, by = 'popA') %>%
    mutate(popA = factor(popA, levels = cluster_pred_order))

  # Data behind the cluster's row of Figure S10 (09_windowed_f3_figure.R).
  write.csv(empirical_results %>%
              dplyr::filter(focal_start == focal_window_for_plot) %>%
              dplyr::transmute(cluster = active_cluster, popA = as.character(popA),
                               ref = as.character(ref), focal_start, focal_z,
                               bg_median_z, bg_q025_z, bg_q975_z,
                               predictor = .data[[cluster_predictor]]),
            sprintf('results/f3/focal_window_f3%s.csv', file_suffix), row.names = FALSE)

  # ---- Focal-window deviation regressed on the cluster's predictor ---------
  #
  # Does focal-window enrichment track the cluster's environmental predictor?
  # Separate models per reference pair.

  focal_central <- empirical_results %>% dplyr::filter(focal_start == focal_window_for_plot)

  cat(sprintf('\n=== Focal-window deviation vs %s, by reference pair (window = %d) ===\n',
              cluster_predictor, focal_window_for_plot))
  focal_lm_fits <- list()
  for (r in c('FL1-CAL1', 'ME1-CAL1', 'ME1-FL1')) {
    dat <- focal_central %>% dplyr::filter(ref == r)
    fit <- lm(reformulate(cluster_predictor, response = 'I(focal_z - bg_median_z)'), data = dat)
    cat('\n---', r, '---\n')
    print(summary(fit))
    focal_lm_fits[[r]] <- fit
  }


  # ===========================================================================
  # MANUSCRIPT VALUES (B): cluster-specific focal-window summary
  # ===========================================================================

  cat('\n\n')
  cat('============================================================\n')
  cat(sprintf('MANUSCRIPT VALUES (B): focal-window f3 [%s]\n', active_cluster))
  cat(sprintf('  Cluster bounds: %s, %.2f-%.2f Mb\n',
              cluster_chr, cluster_start/1e6, cluster_end/1e6))
  cat(sprintf('  Predictor:      %s\n', cluster_predictor))
  cat(sprintf('  Central window: %d\n', focal_window_for_plot))
  cat('============================================================\n')

  # Cluster-window summaries for the central window, per reference pair:
  #   range of focal_z, range of emp_percentile, range of emp_p_z, count of
  #   focal pops falling below the 2.5th background percentile.
  cat(sprintf('\n--- Central window (%d) summary per reference pair ---\n',
              focal_window_for_plot))
  central_summary <- focal_central %>%
    group_by(ref) %>%
    summarize(
      min_focal_z   = min(focal_z),
      median_focal_z= median(focal_z),
      max_focal_z   = max(focal_z),
      min_emp_pct   = min(emp_percentile),
      max_emp_pct   = max(emp_percentile),
      min_emp_p     = min(emp_p_z),
      max_emp_p     = max(emp_p_z),
      n_below_q025  = sum(focal_z < bg_q025_z),
      n_pops        = n(),
      .groups       = 'drop'
    )
  print(central_summary, width = Inf)

  # Full per-population z-scores at the central window for each ref pair.
  cat(sprintf('\n--- Per-population focal_z at central window (%d), sorted ---\n',
              focal_window_for_plot))
  for (r in c('FL1-CAL1', 'ME1-CAL1', 'ME1-FL1')) {
    cat(sprintf('\n  >>> %s <<<\n', r))
    sub <- focal_central %>%
      dplyr::filter(ref == r) %>%
      arrange(focal_z) %>%
      dplyr::select(popA, focal_z, emp_percentile, emp_p_z) %>%
      mutate(
        focal_z        = round(focal_z, 2),
        emp_percentile = round(emp_percentile, 4),
        emp_p_z        = signif(emp_p_z, 3)
      )
    print(sub, n = Inf)
  }

  # Regression of (focal_z - bg_median_z) on the cluster's predictor, per ref pair.
  cat(sprintf('\n--- focal-background deviation vs %s, per reference pair ---\n',
              cluster_predictor))
  for (r in c('FL1-CAL1', 'ME1-CAL1', 'ME1-FL1')) {
    fit <- focal_lm_fits[[r]]
    sm  <- summary(fit)
    co  <- coef(sm)
    pred_row <- cluster_predictor
    cat(sprintf('  %-9s: beta = %7.3f, R^2 = %.3f, adj R^2 = %.3f, p = %.3g\n',
                r,
                co[pred_row, 'Estimate'],
                sm$r.squared,
                sm$adj.r.squared,
                co[pred_row, 'Pr(>|t|)']))
  }

  # Summaries across all overlapping windows (not just the central one), per
  # window x ref pair.
  cat(sprintf('\n--- All cluster windows x reference pairs (min focal_z per cell) ---\n'))
  all_windows_summary <- empirical_results %>%
    group_by(focal_start, ref) %>%
    summarize(
      min_focal_z    = min(focal_z),
      median_focal_z = median(focal_z),
      n_below_q025   = sum(focal_z < bg_q025_z),
      n_pops         = n(),
      .groups        = 'drop'
    )
  print(all_windows_summary, n = Inf)

  cat('\n============================================================\n')
  cat(sprintf('END MANUSCRIPT VALUES (B) -- %s\n', active_cluster))
  cat('============================================================\n\n')

  invisible(empirical_results)
}

for (cl in CLUSTERS) run_cluster(cl)
