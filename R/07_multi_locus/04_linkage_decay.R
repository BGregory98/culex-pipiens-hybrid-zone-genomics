# =============================================================================
# R/07_multi_locus/04_linkage_decay.R
# =============================================================================
# Linkage at the three candidate loci, measured as the pairwise correlation
# (r^2) of allele frequencies across the 24 field pools (reference pools
# excluded; Kapun et al. 2014). Five SNP sets per locus:
#   (A) within-cluster outliers   BayPass candidate SNPs in the cluster
#   (B) within-cluster all-SNPs   PoolSNP-pass SNPs in the cluster, subsampled
#   (C) cluster-flank             outliers x SNPs within 3 Mb outside the cluster
#                                 (block extent and boundary sharpness)
#   (D) cross-locus unlinked      this locus's outliers x the other loci's
#                                 outliers (covariance from shared selection)
#   (E) chr baseline              random non-cluster SNPs on the same chromosome
# Pairs are binned by distance; the binned summary is the table behind
# Figure 4.
#
# Manuscript: Figure 4 (drawn by 05_linkage_decay_figure.R); Results
#             "Cx. quinquefasciatus ancestry and sequencing read depth at a CYP6
#             cluster rise with summer temperature", "Cx. quinquefasciatus
#             ancestry is excluded at a developmental-quiescence locus on the
#             same gradient", "A pre-adapted molestus inversion tracks
#             urbanization independently of temperature"; Supplementary Methods
#             "Linkage disequilibrium analysis"; Supplementary Results "Linkage
#             architecture at the three candidate loci"
#
# Inputs:  results/baypass/v4/baypass_outliers_clustered.csv  (R/03_baypass/01_baypass_outliers.R)
#          data/sample_categories.csv                         (via R/helpers/sample_metadata.R)
#          deposit: data/allele_frequencies_v4.csv            (R/01_filtering/01_allele_frequency_table.R)
# Outputs: results/LD/Table_multi_locus_LD_decay.csv
#
# Runtime: dominated by reading the 1.9 GB allele-frequency table
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")

library(tidyverse)
library(data.table)


# ---- Configuration ---------------------------------------------------------

# Each locus needs: contig, cluster start/end, name used in the table.
# CpL2a is the 500 kb candidate-SNP core here, not the 2.5 Mb block in loci.R.
LD_LOCI <- tribble(
  ~name,    ~contig,  ~start,       ~end,
  'chr1',   CHR1,     CPL1A_START,  CPL1A_END,   # CYP6 cluster
  'chr2',   CHR2,     CPL2A_CORE_START, CPL2A_CORE_END,      # Yan-centred candidate core
  'chr3',   CHR3,     CPL3A_START,  CPL3A_END    # inversion (DELLY breakpoints)
)

# Distance bins (bp). Includes Mb-scale bins to capture chr3's 11.6 Mb
# inversion-scale within-cluster pairs.
DIST_BREAKS <- c(0, 1e3, 5e3, 20e3, 100e3, 250e3, 500e3, 750e3,
                 1e6, 1.5e6, 2e6, 3e6, 5e6, 8e6, 12e6, 20e6, Inf)
DIST_LABELS <- c('<1kb','1-5kb','5-20kb','20-100kb','100-250kb','250-500kb',
                 '500-750kb','750kb-1Mb','1-1.5Mb','1.5-2Mb','2-3Mb',
                 '3-5Mb','5-8Mb','8-12Mb','12-20Mb','>20Mb')

# Subsample sizes (cap matrix size; cor() is O(n^2 * pools))
N_OUTLIERS_MAX       <- 5000   # cap (B) all-SNPs subsample within cluster
N_CHR_BASELINE_SNPS  <- 2000   # (E) random non-cluster baseline per chrom
N_CLUSTER_FLANK_SNPS <- 3000   # (C) flanking SNPs paired with anchors
FLANK_RADIUS         <- 3e6    # (C) how far past the cluster to hunt for the edge

set.seed(20260508)


# ---- Paths -----------------------------------------------------------------

BP_FOCAL_FILE <- 'results/baypass/v4/baypass_outliers_clustered.csv'
need_files(BP_FOCAL_FILE)
AF_FILE       <- deposit_file('data/allele_frequencies_v4.csv')

OUT_DIR   <- 'results/LD'
TABLE_OUT <- file.path(OUT_DIR, 'Table_multi_locus_LD_decay.csv')

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# ---- Load inputs -----------------------------------------------------------

cat('Loading inputs...\n')
bp_focal    <- read.csv(BP_FOCAL_FILE)
af_dt       <- fread(AF_FILE)

missing_in_af <- setdiff(c(field_names, ref_names), colnames(af_dt))
if (length(missing_in_af) > 0) {
  stop('Pools missing from AF table: ', paste(missing_in_af, collapse = ', '))
}

af_wide <- copy(af_dt)
setkey(af_wide, chromosome, position)

cat(sprintf('  bp_focal: %d focal SNPs\n', nrow(bp_focal)))
cat(sprintf('  af_wide:  %d SNPs (wide format)\n', nrow(af_wide)))
cat(sprintf('  field pools: %d, ref pools excluded: %s\n',
            length(field_names), paste(ref_names, collapse = ', ')))


# ---- Helpers ---------------------------------------------------------------

extract_af_mat <- function(positions_df, af_wide, field_names) {
  positions_dt <- as.data.table(positions_df)[, .(chromosome, position)]
  af_sub <- af_wide[positions_dt, on = c('chromosome', 'position'),
                    nomatch = NULL]
  af_sub <- af_sub[order(chromosome, position)]
  
  mat <- as.matrix(af_sub[, ..field_names])
  rownames(mat) <- paste(af_sub$chromosome, af_sub$position, sep = ':')
  
  vars <- apply(mat, 1, var, na.rm = TRUE)
  keep <- which(vars > 0 & !is.na(vars))
  
  list(matrix    = mat[keep, , drop = FALSE],
       positions = af_sub[keep, .(chromosome, position)])
}


compute_pairwise_r2 <- function(af_data, label, locus) {
  mat <- af_data$matrix
  positions <- af_data$positions
  n <- nrow(mat)
  if (n < 2) {
    cat(sprintf('  [%s|%s] too few SNPs (%d) to compute pairs\n',
                locus, label, n))
    return(data.frame())
  }
  cat(sprintf('  [%s|%s] %d SNPs -> %d pairs\n',
              locus, label, n, choose(n, 2)))
  r_mat  <- cor(t(mat), use = 'pairwise.complete.obs')
  r2_mat <- r_mat^2
  ut <- which(upper.tri(r2_mat), arr.ind = TRUE)
  data.frame(
    chromosome_i = positions$chromosome[ut[, 1]],
    position_i   = positions$position[ut[, 1]],
    chromosome_j = positions$chromosome[ut[, 2]],
    position_j   = positions$position[ut[, 2]],
    r2           = r2_mat[ut],
    locus        = locus,
    set          = label
  ) %>%
    mutate(same_chrom = chromosome_i == chromosome_j,
           distance   = ifelse(same_chrom,
                               abs(position_j - position_i),
                               NA_real_))
}


compute_cross_r2 <- function(mat_a, pos_a, mat_b, pos_b, label, locus) {
  if (nrow(mat_a) < 1 || nrow(mat_b) < 1) {
    cat(sprintf('  [%s|%s] empty input\n', locus, label))
    return(data.frame())
  }
  cat(sprintf('  [%s|%s] %d x %d = %d cross-pairs\n',
              locus, label, nrow(mat_a), nrow(mat_b),
              nrow(mat_a) * nrow(mat_b)))
  r_mat <- cor(t(mat_a), t(mat_b), use = 'pairwise.complete.obs')
  expand.grid(i = seq_len(nrow(mat_a)), j = seq_len(nrow(mat_b))) %>%
    mutate(
      chromosome_i = pos_a$chromosome[i],
      position_i   = pos_a$position[i],
      chromosome_j = pos_b$chromosome[j],
      position_j   = pos_b$position[j],
      r2           = c(r_mat^2),
      same_chrom   = chromosome_i == chromosome_j,
      distance     = ifelse(same_chrom,
                            abs(position_j - position_i),
                            NA_real_),
      locus        = locus,
      set          = label
    ) %>%
    dplyr::select(-i, -j)
}


# ---- Per-locus SNP positions (for cross-locus unlinked control) ------------

outlier_positions <- list()
for (l in seq_len(nrow(LD_LOCI))) {
  locus_name <- LD_LOCI$name[l]
  outlier_positions[[locus_name]] <- bp_focal %>%
    filter(chromosome == LD_LOCI$contig[l],
           position  >= LD_LOCI$start[l],
           position  <= LD_LOCI$end[l]) %>%
    dplyr::select(chromosome, position) %>%
    arrange(position)
  cat(sprintf('  %s outliers in cluster: %d\n',
              locus_name, nrow(outlier_positions[[locus_name]])))
}


# =============================================================================
# Per-locus loop
# =============================================================================

all_pairs <- list()

for (l in seq_len(nrow(LD_LOCI))) {
  locus_name <- LD_LOCI$name[l]
  contig     <- LD_LOCI$contig[l]
  start_pos  <- LD_LOCI$start[l]
  end_pos    <- LD_LOCI$end[l]
  
  cat(sprintf('\n========== %s (%s:%g-%g Mb) ==========\n',
              locus_name, contig, start_pos / 1e6, end_pos / 1e6))
  
  # Set A: within-cluster OUTLIERS
  cat(sprintf('\n[%s | A: within-cluster outliers]\n', locus_name))
  outliers <- outlier_positions[[locus_name]]
  if (nrow(outliers) >= 2) {
    outlier_data <- extract_af_mat(outliers, af_wide, field_names)
    all_pairs[[paste0(locus_name, '_A')]] <-
      compute_pairwise_r2(outlier_data, 'A: outliers', locus_name)
  }
  
  # Set B: within-cluster ALL-SNPs subsample
  cat(sprintf('\n[%s | B: within-cluster all-SNPs]\n', locus_name))
  all_in_cluster <- af_wide[chromosome == contig &
                            position  >= start_pos &
                            position  <= end_pos,
                            .(chromosome, position)]
  cat(sprintf('  %d total PoolSNP-pass SNPs in cluster\n', nrow(all_in_cluster)))
  if (nrow(all_in_cluster) > N_OUTLIERS_MAX) {
    all_in_cluster <- all_in_cluster[sample(nrow(all_in_cluster), N_OUTLIERS_MAX)]
  }
  all_in_cluster <- all_in_cluster[order(position)]
  cat(sprintf('  Using %d (subsample cap = %d)\n',
              nrow(all_in_cluster), N_OUTLIERS_MAX))
  all_data <- extract_af_mat(all_in_cluster, af_wide, field_names)
  all_pairs[[paste0(locus_name, '_B')]] <-
    compute_pairwise_r2(all_data, 'B: all-SNPs', locus_name)
  
  # Set E: chr baseline (random non-cluster SNPs on same chrom)
  cat(sprintf('\n[%s | E: chr baseline]\n', locus_name))
  chr_outside <- af_wide[chromosome == contig &
                         (position < start_pos | position > end_pos),
                         .(chromosome, position)]
  baseline_n <- min(N_CHR_BASELINE_SNPS, nrow(chr_outside))
  baseline_sample <- chr_outside[sample(nrow(chr_outside), baseline_n)]
  baseline_sample <- baseline_sample[order(position)]
  baseline_data   <- extract_af_mat(baseline_sample, af_wide, field_names)
  all_pairs[[paste0(locus_name, '_E')]] <-
    compute_pairwise_r2(baseline_data, 'E: chr baseline', locus_name)
  
  # Set C: cluster-flank. Anchors (cluster outliers) x SNPs in a bounded
  # near-cluster flank, so the block-edge zone is densely sampled. (Random
  # genome-wide sampling would put almost no SNPs within a few Mb of the
  # core, which is exactly where the edge is.)
  cat(sprintf('\n[%s | C: cluster-flank]\n', locus_name))
  flank_region <- af_wide[chromosome == contig &
                          ((position >= start_pos - FLANK_RADIUS & position < start_pos) |
                           (position >  end_pos   & position <= end_pos   + FLANK_RADIUS)),
                          .(chromosome, position)]
  flank_n      <- min(N_CLUSTER_FLANK_SNPS, nrow(flank_region))
  flank_sample <- flank_region[sample(nrow(flank_region), flank_n)][order(position)]
  flank_data   <- extract_af_mat(flank_sample, af_wide, field_names)
  if (nrow(outliers) >= 2 && nrow(flank_data$matrix) >= 1) {
    all_pairs[[paste0(locus_name, '_C')]] <- compute_cross_r2(
      outlier_data$matrix,    outlier_data$positions,
      flank_data$matrix,      flank_data$positions,
      'C: cluster-flank', locus_name)
  }
  
  # Set D: cross-locus unlinked (this locus's outliers x other loci's outliers)
  cat(sprintf('\n[%s | D: cross-locus unlinked]\n', locus_name))
  other_loci <- setdiff(LD_LOCI$name, locus_name)
  other_outliers <- bind_rows(outlier_positions[other_loci])
  if (nrow(other_outliers) >= 1 && nrow(outliers) >= 1) {
    other_data <- extract_af_mat(other_outliers, af_wide, field_names)
    all_pairs[[paste0(locus_name, '_D')]] <- compute_cross_r2(
      outlier_data$matrix, outlier_data$positions,
      other_data$matrix,   other_data$positions,
      'D: cross-locus unlinked', locus_name)
  }
}


# =============================================================================
# Distance-binned summary
# =============================================================================

bin_pairs <- function(pairs_df) {
  pairs_df %>%
    mutate(
      distance_bin = cut(distance, DIST_BREAKS,
                         labels = DIST_LABELS,
                         include.lowest = TRUE, right = FALSE),
      distance_bin = ifelse(is.na(distance_bin), 'unlinked',
                            as.character(distance_bin))
    ) %>%
    group_by(locus, set, distance_bin) %>%
    summarize(
      n_pairs   = n(),
      mean_r2   = mean(r2,   na.rm = TRUE),
      median_r2 = median(r2, na.rm = TRUE),
      q25_r2    = quantile(r2, 0.25, na.rm = TRUE),
      q75_r2    = quantile(r2, 0.75, na.rm = TRUE),
      .groups = 'drop'
    )
}

decay_summary <- bind_rows(lapply(all_pairs, bin_pairs))

cat('\n=== DECAY SUMMARY (median r^2 by locus, set, distance) ===\n')
decay_print <- decay_summary %>%
  arrange(locus, set, distance_bin) %>%
  dplyr::select(locus, set, distance_bin, n_pairs, median_r2)
print(as.data.frame(decay_print))

write.csv(decay_summary, TABLE_OUT, row.names = FALSE)


# =============================================================================
# Headline numbers
# =============================================================================

cat('\n=== HEADLINE STATISTICS ===\n')
hl_table <- bind_rows(lapply(all_pairs, function(df) {
  if (nrow(df) == 0) return(NULL)
  data.frame(
    locus     = df$locus[1],
    set       = df$set[1],
    n_pairs   = nrow(df),
    median_r2 = median(df$r2, na.rm = TRUE)
  )
})) %>% arrange(locus, set)
print(hl_table, row.names = FALSE)


# =============================================================================
# Interpretation
# =============================================================================

cat('\n=== INTERPRETATION ===\n')
cat('For each locus, lines (A) and (B) at hundreds of kb scale tell different stories:\n')
cat('  (A) outliers only at high r^2: the BayPass candidate SNPs co-segregate.\n')
cat('  (B) all-SNPs at high r^2:      the whole interval is a haplotype block.\n')
cat('  Both (A) and (B) ~ baseline: independent loci sharing selection, not a block.\n')
cat('  (B) drops below (A) sharply: outliers track a block, but most non-outlier SNPs do not.\n')
cat('\n')
cat('The (B) / (E) ratio at matched distance summarizes haplotype-block evidence:\n')
cat('excess r^2 in the interval over the chromosome-wide null at the same distance,\n')
cat('before any outlier selection.\n')

cat(sprintf('\nTable:  %s\n', TABLE_OUT))
