# =============================================================================
# R/03_baypass/01_baypass_outliers.R
# =============================================================================
# Builds the candidate-SNP table from the BayPass outliers. Reads every SNP that
# was at least a minor outlier for XtXst or for a Bayes factor (classified on
# the cluster), keeps positions in chains of more than ten such outliers with
# gaps <= 10 kb, keeps SNPs that are outliers for both XtXst and BF, and records
# for each SNP the predictor with the highest BF (and the runner-up) and that
# predictor's category. Prints the outlier counts and per-locus predictor
# composition reported in the text.
#
# Manuscript: Methods "Genome--environment association with BayPass" (1,299
#             candidate SNPs); Results "A genome scan identifies three
#             candidate regions with distinct architectures" (654 SNPs at
#             CpL1a, 84.4% bio05; 110 at CpL2a, 81.8% bio05; 64.4% bio05
#             genome-wide);
#             Supplementary Results "Environmental predictor assignment across
#             candidate SNPs". The output table feeds Figure 3 and every
#             locus-level analysis.
#
# Inputs:  deposit: results/baypass/v4/filtered_baypass_data.csv
#                   (cluster/04_baypass/filter_baypass_results.R)
# Outputs: results/baypass/v4/baypass_outliers_clustered.csv
#
# Runtime: about a minute.
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")

library(tidyverse)

FILTERED_FILE <- deposit_file("results/baypass/v4/filtered_baypass_data.csv")
OUT_CLUSTERED <- "results/baypass/v4/baypass_outliers_clustered.csv"


# =============================================================================
# (1) Outlier thresholds
# =============================================================================
#
# Thresholds follow the approach of Leroy et al. (2020, Evol. Appl.
# 13:1171-1188): minor and major outlier classes for XtX and BF, then clusters
# of neighbouring outliers separated by < 10 kb. Here the minor and major XtXst
# thresholds are the 99th and 99.9th percentiles of the POD-calibrated null,
# and the BF thresholds are 10 and 30 dB. The classes (XtX_outlier_type,
# BF_outlier_type: none / minor / major) were assigned on the cluster by
# filter_baypass_results.R; the thresholds are given here only to report them.
# The XtXst values are the POD percentiles written by
# cluster/04_baypass/pod_03_apply_threshold.R (XtX_minor_T.rds, XtX_major_T.rds).

XtX_minor_T <- 71.1249
XtX_major_T <- 106.5979

BF_minor_T  <- 10
BF_major_T  <- 30


# =============================================================================
# (2) Coarse outlier set
# =============================================================================
#
# One row per SNP x predictor for every SNP that is at least a minor outlier
# for XtXst or for BF.

bp_outliers_raw <- read.csv(FILTERED_FILE)

bp_classified <- bp_outliers_raw

cat('\n=== Outlier counts (unique SNPs) ===\n')
cat('Total SNPs with at least one minor outlier type: ',
    length(unique(bp_outliers_raw$MRK)), '\n')
cat('XtXst minor: ', length(unique(bp_outliers_raw$MRK[bp_outliers_raw$XtX_outlier_type == 'minor'])), '\n')
cat('XtXst major: ', length(unique(bp_outliers_raw$MRK[bp_outliers_raw$XtX_outlier_type == 'major'])), '\n')
cat('BF minor: ',    length(unique(bp_outliers_raw$MRK[bp_outliers_raw$BF_outlier_type  == 'minor'])), '\n')
cat('BF major: ',    length(unique(bp_outliers_raw$MRK[bp_outliers_raw$BF_outlier_type  == 'major'])), '\n')


# =============================================================================
# (3) Cumulative genomic coordinates
# =============================================================================
#
# Chromosome lengths used for cumulative offsets (bp):
#   chr1 (ENA_OZ004311_OZ004311.1): 128,522,754
#   chr2 (ENA_OZ004312_OZ004312.1): 213,114,244 (starts at 128,522,754)
#   chr3 (ENA_OZ004313_OZ004313.1):             (starts at 341,636,998)

bp_outliers_raw <- bp_outliers_raw %>%
  mutate(
    cumulative_position = case_when(
      chromosome == CHR1 ~ position,
      chromosome == CHR2 ~ position + 128522754,
      chromosome == CHR3 ~ position + 128522754 + 213114244
    )
  )


# =============================================================================
# (4) Keep outliers in chains of more than ten within 10 kb
# =============================================================================
#
# For each chromosome, outlier positions are chained wherever consecutive
# positions are <= 10 kb apart; chains of more than min_cluster_size positions
# (i.e. >= 11) are kept and all other positions dropped.

cluster_outliers_on_chr <- function(df, chr, max_gap = 10000L, min_cluster_size = 10) {
  # All outlier positions on this chromosome, sorted and unique.
  positions <- sort(unique(df$position[df$chromosome == chr]))
  if (length(positions) < 2L) return(df[FALSE, ])

  # Gaps between consecutive positions (first gap is Inf by convention).
  gaps <- c(Inf, diff(positions))

  # Start a new chain whenever the gap exceeds max_gap.
  cluster_id <- cumsum(gaps > max_gap)
  clusters   <- split(positions, cluster_id)

  clusters_kept      <- clusters[lengths(clusters) > min_cluster_size]
  retained_positions <- unlist(clusters_kept, use.names = FALSE)

  df[df$chromosome == chr & df$position %in% retained_positions, ]
}

chroms <- c(CHR1, CHR2, CHR3)

bp_clustered <- do.call(
  rbind,
  lapply(chroms, cluster_outliers_on_chr, df = bp_outliers_raw)
)

# SNP identifier with ':' as separator (no collision with characters in the
# chromosome names).
bp_clustered$SNP_ID <- factor(
  paste0(bp_clustered$chromosome, ':', bp_clustered$position)
)


# =============================================================================
# (5) Dual outliers and the most strongly associated predictor per SNP
# =============================================================================
#
# The BayPass output has one row per SNP x predictor. Each SNP that is at least
# a minor outlier for both XtXst and BF is reduced to the row of the predictor
# with the highest BF; the second-highest predictor and its BF are kept for
# reference.

# Top-1 and top-2 predictors by BF per SNP.
focal_per_snp <- bp_clustered %>%
  filter(BF_outlier_type != 'none', XtX_outlier_type != 'none') %>%
  group_by(SNP_ID) %>%
  arrange(desc(BF.dB.), .by_group = TRUE) %>%
  summarize(
    predictor_1 = dplyr::first(predictor),
    BF_1        = dplyr::first(BF.dB.),
    predictor_2 = dplyr::nth(predictor, 2),
    BF_2        = dplyr::nth(BF.dB., 2),
    .groups     = 'drop'
  ) %>%
  dplyr::rename(predictor = predictor_1, BF.dB. = BF_1) %>%
  # Recover chromosome and position from SNP_ID.
  separate(
    SNP_ID,
    into    = c('chromosome', 'position'),
    sep     = ':',
    remove  = FALSE,
    convert = TRUE
  ) %>%
  arrange(chromosome, position)

# Join back to the clustered table for XtXst, the outlier classes and the
# cumulative position of the top-predictor row.
bp_focal <- left_join(
  focal_per_snp,
  bp_clustered,
  by = c('SNP_ID', 'chromosome', 'position', 'predictor', 'BF.dB.')
)

bp_focal = bp_focal %>%
  mutate(
    predictor_category = case_when(
      predictor %in% c('bio05', 'bio07') ~ 'temperature',
      predictor %in% c('pd', 'pis', 'tree_canopy') ~ 'urbanization',
      TRUE                                       ~ 'other'
    ))

write.csv(bp_focal, OUT_CLUSTERED, row.names = FALSE)


# =============================================================================
# (6) Predictor composition at the three candidate regions
# =============================================================================
#
# Region bounds (strict inequalities):
#   cluster1  chr1 108.5 - 111.0 Mb  (CpL1a)
#   cluster2  chr2  13.5 -  14.0 Mb  (the 500 kb CpL2a core around Yan)
#   cluster3  chr3 149.40 - 161.02 Mb (CpL3a)

cluster_summaries = bp_focal %>%
  mutate(
    cluster = case_when(
      chromosome == CHR1 & position > CPL1A_START & position < CPL1A_END ~ 'cluster1',
      chromosome == CHR2 & position > CPL2A_CORE_START & position < CPL2A_CORE_END ~ 'cluster2',
      chromosome == CHR3 & position > CPL3A_START & position < CPL3A_END ~ 'cluster3',
      TRUE                                                                ~ 'none'
    )
  ) %>%
  dplyr::count(cluster, predictor, name = 'n_snps') %>%
  arrange(cluster, desc(n_snps))


# =============================================================================
# (7) Values reported in the text
# =============================================================================

cat('\n\n')
cat('============================================================\n')
cat('BayPass outlier identification: values reported in the text\n')
cat('============================================================\n')

cat('\n--- Outlier thresholds ---\n')
cat(sprintf('  XtXst minor (POD 99th percentile):   %.3f\n', XtX_minor_T))
cat(sprintf('  XtXst major (POD 99.9th percentile): %.3f\n', XtX_major_T))
cat(sprintf('  BF minor (dB):                       %.3f\n', BF_minor_T))
cat(sprintf('  BF major (dB):                       %.3f\n', BF_major_T))

cat('\n--- Outlier counts (unique SNPs) ---\n')
xtx_minor_n <- length(unique(bp_classified$MRK[bp_classified$XtX_outlier_type == 'minor']))
xtx_major_n <- length(unique(bp_classified$MRK[bp_classified$XtX_outlier_type == 'major']))
bf_minor_n  <- length(unique(bp_classified$MRK[bp_classified$BF_outlier_type  == 'minor']))
bf_major_n  <- length(unique(bp_classified$MRK[bp_classified$BF_outlier_type  == 'major']))
cat(sprintf('  XtXst outliers: %d total (%d minor + %d major)\n',
            xtx_minor_n + xtx_major_n, xtx_minor_n, xtx_major_n))
cat(sprintf('  BF    outliers: %d total (%d minor + %d major)\n',
            bf_minor_n + bf_major_n, bf_minor_n, bf_major_n))

cat('\n--- Candidate SNPs ---\n')
cat(sprintf('  Total candidate SNPs retained: %d\n', nrow(bp_focal)))

cat('\n--- Genome-wide predictor breakdown (all candidate SNPs) ---\n')
gw_pred <- bp_focal %>% dplyr::count(predictor, name = 'n_snps') %>% arrange(desc(n_snps))
for (i in seq_len(nrow(gw_pred))) {
  cat(sprintf('  %-18s %5d (%5.1f%%)\n',
              as.character(gw_pred$predictor[i]),
              gw_pred$n_snps[i],
              100 * gw_pred$n_snps[i] / nrow(bp_focal)))
}

cat('\n--- Per-region SNP counts and dominant predictors ---\n')
cluster_totals <- cluster_summaries %>%
  filter(cluster != 'none') %>%
  group_by(cluster) %>%
  summarize(total = sum(n_snps), .groups = 'drop')

for (clust in c('cluster1', 'cluster2', 'cluster3')) {
  this_total <- cluster_totals %>% filter(cluster == clust) %>% pull(total)
  if (length(this_total) == 0 || this_total == 0) {
    cat(sprintf('\n  [%s] no SNPs\n', clust))
    next
  }
  cat(sprintf('\n  [%s] total candidate SNPs: %d\n', clust, this_total))
  this_clust <- cluster_summaries %>%
    filter(cluster == clust) %>%
    arrange(desc(n_snps))
  for (i in seq_len(nrow(this_clust))) {
    cat(sprintf('    %-18s %5d (%5.1f%%)\n',
                as.character(this_clust$predictor[i]),
                this_clust$n_snps[i],
                100 * this_clust$n_snps[i] / this_total))
  }
}

cat('\n--- Predictor breakdown across the three regions combined ---\n')
pred_totals <- cluster_summaries %>%
  filter(cluster != 'none') %>%
  group_by(predictor) %>%
  summarize(total = sum(n_snps), .groups = 'drop') %>%
  arrange(desc(total))
grand_total <- sum(pred_totals$total)
cat(sprintf('  (total = %d candidate SNPs in the three regions)\n', grand_total))
for (i in seq_len(nrow(pred_totals))) {
  cat(sprintf('  %-18s %5d (%5.1f%%)\n',
              as.character(pred_totals$predictor[i]),
              pred_totals$total[i],
              100 * pred_totals$total[i] / grand_total))
}

cat('\n============================================================\n\n')
