# =============================================================================
# R/07_multi_locus/06_allele_frequency_clines.R
# =============================================================================
# Allele-frequency clines at the BayPass candidate SNPs of CpL1a (cluster1)
# and CpL2a (cluster2). For each candidate SNP, a binomial GLM of alt/ref
# counts on the cluster's predictor (bio05) across the 24 field pools; the
# slope is then polarized to the cluster's reference pair (CpL1a QUI vs MOL,
# CpL2a QUI vs PIP; R/helpers/cluster_presets.R) and the set restricted to
# SNPs near-fixed between that pair (|freq_A - freq_B| >= fix_threshold).
# Per-pool alternate-allele and total read counts come from the allele-count
# table.
#
# Manuscript: Results "Cx. quinquefasciatus ancestry and sequencing read depth
#             at a CYP6 cluster rise with summer temperature" (544 SNPs, 543 of
#             544); Supplementary Methods "Allele-frequency clines and reference
#             structure at CpL1a and CpL2a"; Supplementary Results "Reference
#             structure at the CpL1a candidate SNPs, and why the QUI-versus-PIP
#             contrast is unavailable" (median slope, P < 0.01 count, PIP median)
#
# Inputs:  results/baypass/v4/baypass_outliers_clustered.csv  (R/03_baypass/01_baypass_outliers.R)
#          data/site_coords_with_spatial_data.csv             (R/00_setup/01_site_predictors.R)
#          data/sample_categories.csv                         (via R/helpers/sample_metadata.R)
#          deposit: data/allele_frequencies_v4.csv            (R/01_filtering/01_allele_frequency_table.R)
#          deposit: data/allele_counts_v4.tsv                 (cluster/04_baypass/vcf_to_allele_counts.sh)
#          deposit: data/vcfs/v4/poolsnp.biallelic.nomiss.vcf.gz  (header only: pool column order)
# Outputs: data/candidate_afs_cluster1.csv, data/candidate_afs_cluster2.csv
#          data/focal_snp_results_cluster1.csv, data/focal_snp_results_cluster2.csv
#            (read by 07_reference_fst_decomposition.R)
#
# Runtime: dominated by reading the 1.9 GB frequency and 2.0 GB count tables
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/cluster_presets.R")
source("R/helpers/sample_metadata.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Clusters whose clines the paper reports
CLUSTERS <- c('cluster1', 'cluster2')

BP_FILE     <- 'results/baypass/v4/baypass_outliers_clustered.csv'
ENV_FILE    <- 'data/site_coords_with_spatial_data.csv'
need_files(BP_FILE, ENV_FILE)
AF_FILE     <- deposit_file('data/allele_frequencies_v4.csv')
COUNTS_FILE <- deposit_file('data/allele_counts_v4.tsv')
VCF_FILE    <- deposit_file('data/vcfs/v4/poolsnp.biallelic.nomiss.vcf.gz')


# -----------------------------------------------------------------------------
# 1. Candidate SNPs per cluster from the BayPass output
# -----------------------------------------------------------------------------

bp_outliers <- read.csv(BP_FILE)

candidate_snps_by_cluster <- lapply(setNames(CLUSTERS, CLUSTERS), function(cl) {
  cfg <- cluster_presets[[cl]]
  bp_outliers %>%
    filter(chromosome == cfg$chr,
           position   >  cfg$start,
           position   <  cfg$end) %>%
    dplyr::select(chromosome, position)
})
all_candidate_snps <- bind_rows(candidate_snps_by_cluster)


# -----------------------------------------------------------------------------
# 2. Allele frequencies and depths at the candidate SNPs (read once)
# -----------------------------------------------------------------------------

recode_chr <- function(x) recode(x,
  'ENA_OZ004311_OZ004311.1' = 'Chr1',
  'ENA_OZ004312_OZ004312.1' = 'Chr2',
  'ENA_OZ004313_OZ004313.1' = 'Chr3'
)

# fread rather than read.csv: the table is 1.9 GB and read.csv exceeds R's
# vector limit before the semi_join cuts it to the candidate SNPs.
af_candidates <- as.data.frame(data.table::fread(AF_FILE, showProgress = FALSE)) %>%
  semi_join(all_candidate_snps, by = c('chromosome', 'position'))

# Sample names in VCF column order, read from the #CHROM header line; the
# count columns of the allele-count table follow the same order.
samples <- local({
  con <- gzfile(VCF_FILE, 'rt'); on.exit(close(con))
  repeat {
    ln <- readLines(con, n = 1)
    if (!length(ln)) stop('no #CHROM header line found in ', VCF_FILE)
    if (startsWith(ln, '#CHROM')) break
  }
  strsplit(ln, '\t', fixed = TRUE)[[1]][-(1:9)]
})

# Read counts table and parse "RD:AD" format into total read depth
counts_wide <- as.data.frame(data.table::fread(COUNTS_FILE,
                                               header = FALSE,
                                               showProgress = FALSE)) %>%
  dplyr::rename_with(~ paste0('X', seq_along(.))) %>%
  dplyr::rename(chromosome = X1, position = X2, ref = X3, alt = X4) %>%
  mutate(chromosome = recode_chr(chromosome),
         snp_id = paste0(chromosome, '-', position))

# Restrict to the candidate SNPs before parsing. counts_wide is only ever joined
# to the candidate set, so this is identical downstream and avoids parsing
# 11.6M rows.
counts_wide <- counts_wide %>%
  dplyr::filter(snp_id %in% paste0(recode_chr(af_candidates$chromosome), '-',
                                   af_candidates$position))

# Parse RD:AD columns into total depth (counts_wide) and alternate-allele
# depth (alt_wide) for each sample
alt_wide <- counts_wide[, 'snp_id', drop = FALSE]
for (i in seq_along(samples)) {
  col   <- counts_wide[[4 + i]]
  parts <- str_split_fixed(col, ':', 2)
  counts_wide[[samples[i]]] <- as.integer(parts[, 1]) + as.integer(parts[, 2])
  alt_wide[[samples[i]]]    <- as.integer(parts[, 2])
}

site_env <- read.csv(ENV_FILE) %>%
  dplyr::rename(sample = site) %>%
  dplyr::select(sample, bio05, bio07, bio18, elevation, pd, pis, tree_canopy, wtd)


# =============================================================================
# Per-cluster analysis
# =============================================================================

run_cluster <- function(active_cluster) {

  cfg                <- cluster_presets[[active_cluster]]
  cluster_chr        <- cfg$chr
  cluster_start      <- cfg$start
  cluster_end        <- cfg$end
  polarize_ref_A     <- cfg$polarize_A
  polarize_ref_B     <- cfg$polarize_B
  fix_threshold      <- cfg$fix_threshold
  cluster_predictor  <- cfg$predictor
  cluster_pred_label <- cfg$predictor_label
  file_suffix        <- paste0('_', active_cluster)

  cat(sprintf('\nActive cluster: %s  (%s, %.2f-%.2f Mb)\n',
              active_cluster, cluster_chr,
              cluster_start / 1e6, cluster_end / 1e6))
  cat(sprintf('Polarizing %s (target, y-axis) vs. %s (contrast)\n',
              polarize_ref_A, polarize_ref_B))
  cat(sprintf('Predictor:  %s  (%s)\n',
              cluster_predictor, cluster_pred_label))
  cat(sprintf('Fixation filter: |freq_A - freq_B| >= %.2f\n\n', fix_threshold))

  candidate_snps <- candidate_snps_by_cluster[[active_cluster]]
  cat('Candidate SNPs in cluster:', nrow(candidate_snps), '\n')


  # ---------------------------------------------------------------------------
  # 3. Allele frequency table for the cluster's candidate SNPs
  # ---------------------------------------------------------------------------

  af_long <- af_candidates %>%
    semi_join(candidate_snps, by = c('chromosome', 'position')) %>%
    dplyr::select(-any_of(c('region_type', 'is_coding'))) %>%
    mutate(chromosome = recode_chr(chromosome)) %>%
    pivot_longer(
      cols      = -c(chromosome, position, cumulative_position),
      names_to  = 'sample',
      values_to = 'frequency'
    ) %>%
    filter(!sample %in% c('Northfield', 'Evanston', 'FL2', 'FL3', 'ME2', 'ME3')) %>%
    mutate(snp_id = paste0(chromosome, '-', position))


  # ---------------------------------------------------------------------------
  # 4. Separate reference and field population allele frequencies
  # ---------------------------------------------------------------------------

  ref_afs <- af_long %>%
    filter(sample %in% ref_names) %>%
    dplyr::select(snp_id, sample, frequency) %>%
    pivot_wider(names_from = sample, values_from = frequency) %>%
    dplyr::rename(ME1_freq = ME1, FL1_freq = FL1, CAL1_freq = CAL1)

  candidate_afs <- af_long %>%
    filter(!sample %in% ref_names) %>%
    left_join(site_env, by = 'sample')


  # ---------------------------------------------------------------------------
  # 5. Ref/alt counts per field sample
  # ---------------------------------------------------------------------------

  counts_cl <- counts_wide %>% dplyr::filter(snp_id %in% unique(af_long$snp_id))

  # Extract ref/alt alleles for reference populations
  ref_alleles <- counts_cl %>%
    dplyr::select(snp_id, ref, alt, all_of(ref_names))

  ref_afs <- ref_afs %>%
    left_join(ref_alleles, by = 'snp_id') %>%
    dplyr::select(snp_id, ref, alt, ME1_freq, FL1_freq, CAL1_freq)

  # Compute alt and ref counts for field populations
  candidate_afs <- candidate_afs %>%
    left_join(
      counts_cl %>%
        dplyr::select(snp_id, all_of(field_names)) %>%
        pivot_longer(-snp_id, names_to = 'sample', values_to = 'count'),
      by = c('snp_id', 'sample')
    ) %>%
    left_join(
      alt_wide %>%
        dplyr::filter(snp_id %in% unique(af_long$snp_id)) %>%
        dplyr::select(snp_id, all_of(field_names)) %>%
        pivot_longer(-snp_id, names_to = 'sample', values_to = 'alt_ct'),
      by = c('snp_id', 'sample')
    ) %>%
    mutate(ref_ct = count - alt_ct) %>%
    dplyr::select(snp_id, chromosome, position, cumulative_position, sample,
           count, alt_ct, ref_ct, frequency,
           bio05, bio07, bio18, elevation, pd, pis, tree_canopy, wtd)

  write.csv(candidate_afs,
            sprintf('data/candidate_afs%s.csv', file_suffix),
            row.names = FALSE)


  # ---------------------------------------------------------------------------
  # 6. Binomial GLMs: regress alt allele frequency on the cluster's predictor
  # ---------------------------------------------------------------------------
  #
  # GLM formula is built from cluster_predictor via reformulate().

  snp_ids <- unique(na.omit(candidate_afs$snp_id))

  glm_formula <- reformulate(cluster_predictor, response = 'cbind(alt_ct, ref_ct)')

  glm_results <- map_dfr(snp_ids, function(id) {
    sub   <- candidate_afs[candidate_afs$snp_id == id, ]
    fit   <- glm(glm_formula, family = binomial, data = sub)
    coefs <- summary(fit)$coefficients
    # Row 2 = the single predictor; keep position-based indexing to avoid
    # name-matching issues when cluster_predictor is a column name.
    data.frame(snp_id = id, slope = coefs[2, 1], p_val = coefs[2, 4])
  })


  # ---------------------------------------------------------------------------
  # 7. Polarize SNPs relative to the chosen reference pair + fixation filter
  # ---------------------------------------------------------------------------
  #
  # For each SNP, identify which of polarize_ref_A / polarize_ref_B has the
  # higher alt-allele frequency. Then:
  #   - `higher_af`   records which reference is "alt-like" for that SNP
  #   - `ref_af_diff` is signed as (freq_A - freq_B)
  #   - `slope`       is flipped so positive always means "increasing
  #                   polarize_ref_A-like allele with the predictor"
  #
  # After polarization, restrict to SNPs that are effectively fixed between the
  # two polarization references (|ref_af_diff| >= fix_threshold). At SNPs where
  # polarize_ref_A and polarize_ref_B carry similar frequencies, the
  # polarization assignment is effectively arbitrary and adds noise to the
  # slope distribution. The threshold is set per cluster in the preset (0.95 at
  # CpL1a; 0.70 at CpL2a, where the references are less differentiated).

  freq_col_A <- paste0(polarize_ref_A, '_freq')
  freq_col_B <- paste0(polarize_ref_B, '_freq')

  results_unfiltered <- glm_results %>%
    left_join(ref_afs, by = 'snp_id') %>%
    mutate(
      freq_A      = .data[[freq_col_A]],
      freq_B      = .data[[freq_col_B]],
      higher_af   = case_when(
        freq_A > freq_B ~ polarize_ref_A,
        freq_B > freq_A ~ polarize_ref_B,
        TRUE            ~ 'equal'
      ),
      ref_af_diff = freq_A - freq_B,
      # Flip slope when alt is polarize_ref_B-like, so that positive slope
      # consistently means "increasing polarize_ref_A-like allele".
      slope       = ifelse(higher_af == polarize_ref_B, -slope, slope)
    )

  # Fixation filter.
  results <- results_unfiltered %>%
    filter(abs(ref_af_diff) >= fix_threshold)

  n_pre  <- nrow(results_unfiltered)
  n_post <- nrow(results)

  cat(sprintf(
    '\nFixation filter (|ref_af_diff| >= %.2f, %s vs %s):\n',
    fix_threshold, polarize_ref_A, polarize_ref_B
  ))
  cat(sprintf('  SNPs before filter: %d\n',  n_pre))
  cat(sprintf('  SNPs after  filter: %d  (%.1f%% retained)\n',
              n_post, 100 * n_post / n_pre))

  if (n_post < 50) {
    warning(sprintf(
      'Only %d SNPs retained after fixation filter at threshold %.2f. ',
      n_post, fix_threshold
    ), 'Consider lowering fix_threshold in the cluster preset.')
  }

  write.csv(results,
            sprintf('data/focal_snp_results%s.csv', file_suffix),
            row.names = FALSE)

  cat(sprintf('\nFocal SNPs retained (polarized %s vs %s): %d\n',
              polarize_ref_A, polarize_ref_B, nrow(results)))
  cat('Slope summary (filtered) -- median:', round(median(results$slope), 3),
      ' min:', round(min(results$slope), 3),
      ' max:', round(max(results$slope), 3), '\n')
  cat('Proportion of SNPs with P < 0.001:',
      round(mean(results$p_val < 0.001), 3), '\n')
  cat('Proportion with slope > 0 (i.e., increasing ',
      polarize_ref_A, '-like allele): ',
      round(mean(results$slope > 0), 3), '\n\n', sep = '')


  # ---------------------------------------------------------------------------
  # 8. Reference medians on the polarized scale
  # ---------------------------------------------------------------------------
  #
  # Median allele frequency of each reference, polarized to polarize_ref_A.

  ref_medians <- results %>%
    mutate(
      FL1_polar  = ifelse(higher_af == polarize_ref_B, 1 - FL1_freq,  FL1_freq),
      ME1_polar  = ifelse(higher_af == polarize_ref_B, 1 - ME1_freq,  ME1_freq),
      CAL1_polar = ifelse(higher_af == polarize_ref_B, 1 - CAL1_freq, CAL1_freq)
    ) %>%
    summarise(
      FL1_med  = median(FL1_polar,  na.rm = TRUE),
      ME1_med  = median(ME1_polar,  na.rm = TRUE),
      CAL1_med = median(CAL1_polar, na.rm = TRUE)
    )


  # ===========================================================================
  # MANUSCRIPT VALUES: candidate-AF cline summary for the active cluster
  # ===========================================================================

  cat('\n\n')
  cat('============================================================\n')
  cat(sprintf('MANUSCRIPT VALUES: candidate-AF clines [%s]\n', active_cluster))
  cat(sprintf('  Cluster bounds: %s, %.2f-%.2f Mb\n',
              cluster_chr, cluster_start/1e6, cluster_end/1e6))
  cat(sprintf('  Polarization:   target=%s, contrast=%s\n',
              polarize_ref_A, polarize_ref_B))
  cat(sprintf('  Predictor:      %s\n', cluster_predictor))
  cat(sprintf('  Fixation filter: |freq_A - freq_B| >= %.2f\n',
              fix_threshold))
  cat('============================================================\n')

  # --- Sample sizes before/after the fixation filter ---
  cat('\n--- SNP counts ---\n')
  cat(sprintf('  Candidate SNPs in cluster (from BayPass):            %d\n',
              nrow(candidate_snps)))
  cat(sprintf('  SNPs with fitted binomial GLM:                       %d\n',
              nrow(glm_results)))
  cat(sprintf('  SNPs before fixation filter:                         %d\n', n_pre))
  cat(sprintf('  SNPs after  fixation filter:                         %d  (%.1f%%)\n',
              n_post, 100 * n_post / n_pre))

  # --- Slope distribution on the filtered set ---
  cat('\n--- Slope distribution on filtered SNPs (log-odds, polarized) ---\n')
  cat(sprintf('  median:      %8.3f\n', median(results$slope)))
  cat(sprintf('  min:         %8.3f\n', min(results$slope)))
  cat(sprintf('  max:         %8.3f\n', max(results$slope)))
  cat(sprintf('  mean:        %8.3f\n', mean(results$slope)))
  cat(sprintf('  sd:          %8.3f\n', sd(results$slope)))

  # --- Directionality and significance proportions ---
  cat('\n--- Proportions across filtered SNPs ---\n')
  cat(sprintf('  P < 0.001:               %5.1f%% (%d / %d)\n',
              100 * mean(results$p_val < 0.001),
              sum(results$p_val  < 0.001), nrow(results)))
  cat(sprintf('  P < 0.01:                %5.1f%% (%d / %d)\n',
              100 * mean(results$p_val < 0.01),
              sum(results$p_val  < 0.01), nrow(results)))
  cat(sprintf('  slope > 0 (increasing %s-like): %5.1f%% (%d / %d)\n',
              polarize_ref_A,
              100 * mean(results$slope > 0),
              sum(results$slope > 0), nrow(results)))
  cat(sprintf('  slope > 0 AND P < 0.001:         %5.1f%% (%d / %d)\n',
              100 * mean(results$slope > 0 & results$p_val < 0.001),
              sum(results$slope > 0 & results$p_val < 0.001), nrow(results)))

  # --- Reference medians on the polarized scale ---
  cat('\n--- Reference median allele frequencies (polarized to target) ---\n')
  cat(sprintf('  %s (target):   %.3f\n', polarize_ref_A,
              ref_medians[[paste0(polarize_ref_A, '_med')]]))
  cat(sprintf('  %s (contrast): %.3f\n', polarize_ref_B,
              ref_medians[[paste0(polarize_ref_B, '_med')]]))
  third_ref <- setdiff(c('FL1', 'ME1', 'CAL1'), c(polarize_ref_A, polarize_ref_B))
  cat(sprintf('  %s (third):    %.3f\n', third_ref,
              ref_medians[[paste0(third_ref, '_med')]]))

  cat('\n============================================================\n')
  cat(sprintf('END MANUSCRIPT VALUES -- %s\n', active_cluster))
  cat('============================================================\n\n')

  invisible(results)
}

for (cl in CLUSTERS) run_cluster(cl)
