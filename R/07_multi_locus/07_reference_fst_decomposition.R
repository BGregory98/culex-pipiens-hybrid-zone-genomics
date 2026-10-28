# =============================================================================
# R/07_multi_locus/07_reference_fst_decomposition.R
# =============================================================================
# Reference structure at the near-fixed candidate SNPs of CpL1a (cluster1) and
# CpL2a (cluster2): the BayPass candidate SNPs in the cluster that pass the
# fixation filter of 06_allele_frequency_clines.R. For each pair of reference
# pools, FST (poolfstat computeFST, Anova method) restricted to those SNPs,
# with block-jackknife standard errors and 95% CIs from blocks of about 10% of
# the SNPs; plus the variance explained by PC1/PC2 of a PCA of all 27 pools at
# those SNPs.
#
# Manuscript: Supplementary Methods "Allele-frequency clines and reference
#             structure at CpL1a and CpL2a"; Supplementary Results "Reference
#             structure at the CpL1a candidate SNPs, and why the QUI-versus-PIP
#             contrast is unavailable" (CpL1a pairwise FST); Supplementary
#             Methods "Donor-private ancestry panels and the per-population
#             donor-ancestry test" (CpL2a FST(PIP, MOL) close to 0)
#
# Inputs:  results/baypass/v4/baypass_outliers_clustered.csv  (R/03_baypass/01_baypass_outliers.R)
#          data/focal_snp_results_cluster1.csv, data/focal_snp_results_cluster2.csv
#                                                             (06_allele_frequency_clines.R)
#          deposit: data/vcfs/v4/poolsnp.biallelic.nomiss.vcf.gz
# Outputs: console only (pairwise FST and CIs per reference pair, PCA variance)
#
# Runtime: dominated by vcf2pooldata reading the whole VCF
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/cluster_presets.R")

suppressPackageStartupMessages({
  library(poolfstat)
  library(tidyverse)
  library(VariantAnnotation)
})

# Clusters whose reference structure the paper reports
CLUSTERS <- c('cluster1', 'cluster2')

BP_FILE <- 'results/baypass/v4/baypass_outliers_clustered.csv'
FOCAL_SNP_FILES <- setNames(sprintf('data/focal_snp_results_%s.csv', CLUSTERS), CLUSTERS)
need_files(BP_FILE, FOCAL_SNP_FILES)
VCF_FILE <- deposit_file('data/vcfs/v4/poolsnp.biallelic.nomiss.vcf.gz')


# ---- 1. Pool data from the VCF (read once) ----------------------------------

# Pool names in VCF column order, and a lookup of each pool's column index.
names       <- samples(scanVcfHeader(VCF_FILE))
names_index <- data.frame(names, index = seq_along(names))

pooldata <- vcf2pooldata(
  vcf.file  = VCF_FILE,
  poolsizes = rep(80, 27),
  poolnames = names
)

# ---- 2. SNP metadata -------------------------------------------------------

snp_info <- as.data.frame(pooldata@snp.info)
colnames(snp_info) <- tolower(colnames(snp_info))
snp_info$snp_index <- seq_len(nrow(snp_info))
snp_info$chromosome <- gsub("\\|", "_", snp_info$chromosome)

bp_outliers <- read.csv(BP_FILE)

# ---- Population indices for the three references ---------------------------

idx_FL1  <- names_index$index[names_index$names == 'FL1']
idx_ME1  <- names_index$index[names_index$names == 'ME1']
idx_CAL1 <- names_index$index[names_index$names == 'CAL1']

# ---- Reference pairs to compare --------------------------------------------

ref_pairs <- list(
  list(name_A = 'FL1',  name_B = 'CAL1', idx_A = idx_FL1,  idx_B = idx_CAL1),
  list(name_A = 'ME1',  name_B = 'CAL1', idx_A = idx_ME1,  idx_B = idx_CAL1),
  list(name_A = 'ME1',  name_B = 'FL1',  idx_A = idx_ME1,  idx_B = idx_FL1)
)


# =============================================================================
# Per-cluster analysis
# =============================================================================

run_cluster <- function(active_cluster) {

  cfg            <- cluster_presets[[active_cluster]]
  cluster_chr    <- cfg$chr
  cluster_start  <- cfg$start
  cluster_end    <- cfg$end

  cat(sprintf('\nActive cluster: %s  (%s, %.2f-%.2f Mb)\n\n',
              active_cluster, cluster_chr,
              cluster_start / 1e6, cluster_end / 1e6))

  # ---- 3. Identify focal cluster SNPs from BayPass outliers ----------------
  #
  # Two-step filter:
  #   (a) BayPass outlier SNPs falling within the cluster's genomic bounds.
  #   (b) Restrict (a) to the near-fixation subset written by
  #       06_allele_frequency_clines.R (|freq_A - freq_B| >= fix_threshold
  #       between polarize_A and polarize_B).
  # This is the SNP set behind the polarized allele-frequency cline, so the
  # FST decomposition describes the same loci as the cline.

  # (a) BayPass outliers within cluster bounds.
  cluster_snps <- bp_outliers %>%
    filter(chromosome == cluster_chr,
           position   >  cluster_start,
           position   <  cluster_end) %>%
    dplyr::select(chromosome, position)

  # (b) Near-fixation filter. The focal_snp_results file is per cluster, so all
  # rows are on cluster_chr and only `position` is needed for the intersect.
  fixation_filtered_positions <- read.csv(FOCAL_SNP_FILES[[active_cluster]]) %>%
    mutate(position = sapply(strsplit(as.character(snp_id), "-"), `[`, 2)) %>%
    pull(position)

  candidate_snps <- cluster_snps %>%
    filter(position %in% fixation_filtered_positions)

  focal_snp_indices <- snp_info %>%
    inner_join(candidate_snps, by = c('chromosome', 'position')) %>%
    pull(snp_index)

  cat('Cluster SNPs (BayPass outliers within bounds):', nrow(cluster_snps),      '\n')
  cat('Near-fixation-filtered SNPs:                  ', nrow(candidate_snps),    '\n')
  cat('Focal SNPs matched in pooldata:               ', length(focal_snp_indices), '\n\n')


  # ---- 4. Compute FST for each reference pair at focal SNPs ----------------

  fst_results <- map_dfr(ref_pairs, function(pair) {

    pd_sub <- pooldata.subset(
      pooldata,
      pool.index = c(pair$idx_A, pair$idx_B),
      snp.index  = focal_snp_indices,
      verbose    = FALSE
    )

    fst <- computeFST(
      pd_sub,
      method = 'Anova',
      nsnp.per.bjack.block = round(length(focal_snp_indices) / 10, 0),
      verbose = FALSE
    )

    fst_df <- as.data.frame(fst$Fst)

    data.frame(
      pair        = paste0(pair$name_A, '_vs_', pair$name_B),
      popA        = pair$name_A,
      popB        = pair$name_B,
      estimate    = fst_df[1, 1],
      bjack_mean  = fst_df[2, 1],
      bjack_se    = fst_df[3, 1],
      CI95inf     = fst_df[4, 1],
      CI95sup     = fst_df[5, 1],
      n_snps      = length(focal_snp_indices)
    )
  })

  print(fst_results)


  # ---- 5. PCA of focal-SNP allele frequencies, all populations -------------
  #
  # The focal SNPs are near-fixed between polarize_A and polarize_B, so PC1 is
  # dominated by that contrast and PC2 captures variation orthogonal to it.

  # (a) Allele frequencies at focal SNPs (rows = SNPs, cols = pools).
  ref_counts <- pooldata@refallele.readcount[focal_snp_indices, , drop = FALSE]
  tot_counts <- pooldata@readcoverage     [focal_snp_indices, , drop = FALSE]
  af_focal   <- ref_counts / tot_counts
  colnames(af_focal) <- names

  # poolsnp.biallelic.nomiss should have no zero-coverage cells; any that slip
  # through are replaced with the column mean so prcomp does not fail.
  if (any(!is.finite(af_focal))) {
    for (j in seq_len(ncol(af_focal))) {
      bad <- !is.finite(af_focal[, j])
      if (any(bad)) af_focal[bad, j] <- mean(af_focal[!bad, j], na.rm = TRUE)
    }
  }

  # (b) PCA. Populations as rows, SNPs as variables. Mean-center but do not
  # scale -- allele frequencies are already on a common [0,1] scale, and
  # scaling each SNP would over-weight near-monomorphic SNPs (which is most
  # of the focal set after the fixation filter).
  pca     <- prcomp(t(af_focal), scale. = FALSE)
  var_pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)


  # ===========================================================================
  # MANUSCRIPT VALUES: FST decomposition + focal-SNP PCA for the active cluster
  # ===========================================================================

  cat('\n\n')
  cat('============================================================\n')
  cat(sprintf('MANUSCRIPT VALUES: FST decomposition [%s]\n', active_cluster))
  cat(sprintf('  Cluster bounds:        %s, %.2f-%.2f Mb\n',
              cluster_chr, cluster_start/1e6, cluster_end/1e6))
  cat(sprintf('  Cluster SNPs (pre-fixation filter):  %d\n', nrow(cluster_snps)))
  cat(sprintf('  Focal SNPs (post-fixation filter):   %d\n', length(focal_snp_indices)))
  cat('============================================================\n\n')

  # Pairwise FST table with point estimate, block-jackknife mean, SE and 95% CI.
  cat('Pairwise FST among reference populations at focal SNPs:\n')
  cat(sprintf('%-14s %10s %10s %10s %10s %10s\n',
              'pair', 'estimate', 'bjack', 'SE', 'CI95_lo', 'CI95_hi'))
  for (i in seq_len(nrow(fst_results))) {
    r <- fst_results[i, ]
    cat(sprintf('%-14s %10.3f %10.3f %10.4f %10.3f %10.3f\n',
                r$pair, r$estimate, r$bjack_mean, r$bjack_se,
                r$CI95inf, r$CI95sup))
  }

  top_pair <- fst_results[which.max(fst_results$bjack_mean), ]
  cat(sprintf('\nMost differentiated pair: %s  (FST = %.3f [%.3f-%.3f])\n',
              top_pair$pair, top_pair$bjack_mean,
              top_pair$CI95inf, top_pair$CI95sup))

  # Focal-SNP PCA summary.
  cat('\nFocal-SNP PCA (all pools):\n')
  cat(sprintf('  PC1 variance explained:  %.1f%%\n', var_pct[1]))
  cat(sprintf('  PC2 variance explained:  %.1f%%\n', var_pct[2]))
  cat(sprintf('  PC1 + PC2:               %.1f%%\n', var_pct[1] + var_pct[2]))
  cat(sprintf('  Polarization axis:       %s -- %s\n',
              cfg$polarize_A, cfg$polarize_B))

  cat('\n============================================================\n')
  cat(sprintf('END MANUSCRIPT VALUES -- %s\n', active_cluster))
  cat('============================================================\n\n')

  invisible(fst_results)
}

for (cl in CLUSTERS) run_cluster(cl)
