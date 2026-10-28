# =============================================================================
# cluster/07_multilocus/01_genome_wide_null.R
# =============================================================================
# Empirical null distribution of the locus-level statistics: for each
# candidate locus, the same five statistics computed at random autosomal
# regions of the locus's width. Submitted through genome_wide_null.slurm
# (10,000 regions per locus, seed 42).
#
# Step: 9, genome-wide null (cluster/README.md)
# Manuscript: Methods "Locus-level diagnostics against a genome-wide null";
#             Supplementary Methods "Genome-wide null distribution for
#             locus-level diagnostic statistics"
#
# Inputs:  --vcf PoolSNP_out_v4/poolsnp.biallelic.nomiss.vcf.gz (step 5)
#          --env_file data/metadata/site_coords_with_spatial_data.csv
# Outputs: --out_rds genome_wide_null.rds (repo
#          results/multilocus/genome_wide_null.rds; read by
#          R/07_multi_locus/01_per_pool_ancestry.R and
#          03_haplotype_block_diagnostics.R)
#
# Runtime: ~5 h at 10,000 regions per locus (single core, 64 GB).
#
# CONCEPT
# -------
# For each candidate locus (chr1, chr2, chr3), the focal-locus test asks "is
# this statistic exceptional vs random regions of the same width?" Width is
# matched per locus to the haplotype-block scale:
#   chr1 random regions = 2.5 Mb (haplotype block 108.5-111.0 Mb)
#   chr2 random regions = 2.5 Mb (haplotype block 12.5-15.0 Mb)
#   chr3 random regions = 11.62 Mb (inversion 149,402,438-161,021,387, DELLY)
# Each locus thus has its own null distribution. The chr3 width is large;
# non-overlapping autosomal 12-Mb windows number ~250, so random chr3 regions
# overlap, which is fine for an empirical percentile null.
#
# Random regions are drawn from autosomes EXCLUDING the focal regions themselves
# (chr1 108.5-111 Mb, chr2 12.5-15 Mb, chr3 inversion 149.4-161.0 Mb), the same
# exclusion that defines the genome-wide reference panel.
#
# STATISTICS COMPUTED PER RANDOM REGION
# -------------------------------------
# For each random region, per locus configuration (= donor identity + width),
# we compute and store:
#
#   (1) Delta_marginal: median over pools of [alpha_D(region) - alpha_D(genome)]
#   (2) Delta_slope:    slope of [alpha_D(region) - alpha_D(genome)] ~ env
#   (3) PC1_pct_var:    PC1 % variance explained by per-pool AF matrix
#                       restricted to region SNPs
#   (4) PC1_FST_cor:    |cor(PC1 score, mean per-pool FST to donor reference)|
#   (5) FST_donor_var:  variance over pools of mean FST(pool, donor) at region
#
# The chr2 locus is special: we compute (1) and (2) for BOTH ME1-private and
# non-FL1-deficit panels. So the chr2 null carries a "dual" donor
# configuration.
#
# OUTPUT: genome_wide_null.rds  (a named list)
# --------------------------------------------
#   $config            : run parameters (n_random_regions, seed, panel cutoffs)
#   $panels            : list with donor-private SNP indices for FL1, ME1,
#                        CAL1, plus non_FL1_deficit (chr2 dual-panel)
#   $genome_panel_alpha: for each donor, vector length 24 (field pools) of
#                        alpha_D,i computed on the genome-wide panel
#                        (= panel SNPs NOT in any focal broad region)
#   $null$<locus>      : data.frame with one row per random region, columns
#                        for the 5 statistics above, plus n_snps and metadata
#                        (For chr2: columns are present for both panels.)
#
# DOWNSTREAM USE
# --------------
# R/07_multi_locus/01_per_pool_ancestry.R reads $genome_panel_alpha for its
# baseline alpha_genome values, and reads $null$<locus> to compute
# permutation p-values for the focal Delta statistics.
#
# R/07_multi_locus/03_haplotype_block_diagnostics.R reads $null$<locus>
# columns PC1_pct_var, PC1_FST_cor, FST_donor_var to compute permutation
# p-values for the three block-structure statistics.
#
# SCALING
# -------
# The script default is --n_regions 1000; the published run passes 10000
# from the SLURM wrapper. Stage 1's vcfR read of the full VCF (~11.6M SNPs x
# 27 pools) takes ~30-60 minutes. Stage 5 (the random-region loop) dominates
# and scales linearly with n_random_regions; chr3 is the longest since its
# 11.6-Mb regions hold ~5x more SNPs than 2.5-Mb ones.
#
# =============================================================================

# ---- Libraries --------------------------------------------------------------

suppressPackageStartupMessages({
  library(vcfR)
  library(data.table)
  library(optparse)
})


# ---- CLI options (defaults; the SLURM wrapper overrides) -------------------

option_list <- list(
  make_option('--base_dir', type = 'character', default = '.',
              help = 'Project root [default %default]'),
  make_option('--vcf', type = 'character',
              default = 'SequencingData/data/baypass_v4/poolsnp.biallelic.nomiss.vcf.gz',
              help = 'Path to v4 VCF [default %default]'),
  make_option('--env_file', type = 'character',
              default = 'SequencingData/data/site_coords_with_spatial_data.csv',
              help = 'Env covariates CSV [default %default]'),
  make_option('--out_rds', type = 'character',
              default = 'SequencingData/results/multilocus/genome_wide_null.rds',
              help = 'Output RDS path [default %default]'),
  make_option('--n_regions', type = 'integer', default = 1000L,
              help = 'Number of random regions per locus [default %default]'),
  make_option('--seed', type = 'integer', default = 42L,
              help = 'Random seed [default %default]'),
  make_option('--min_panel_snps', type = 'integer', default = 5L,
              help = 'Min panel SNPs per region to retain [default %default]')
)
opt <- parse_args(OptionParser(option_list = option_list))


# ---- Configuration ----------------------------------------------------------

BASE_DIR <- opt$base_dir
setwd(BASE_DIR)

VCF_FILE <- opt$vcf
ENV_FILE <- opt$env_file
OUT_RDS  <- opt$out_rds

dir.create(dirname(OUT_RDS), recursive = TRUE, showWarnings = FALSE)

# Run parameters
n_random_regions <- opt$n_regions
random_seed      <- opt$seed

# Donor-private SNP panel cutoffs
PANEL_P_DONOR <- 0.95
PANEL_P_OTHER <- 0.05

# Minimum depth for a pool AF to count
MIN_DEPTH <- 10

# Minimum SNPs in a random region's panel for the region to contribute to the
# null. Regions with too few panel SNPs produce noisy alpha estimates; we drop
# them rather than letting them inflate null variance.
MIN_PANEL_SNPS_REGION <- opt$min_panel_snps

# Reference pools (excluded from per-pool tests; used to define donor panels)
REF_POOLS <- c('FL1', 'ME1', 'CAL1')

cat('Configuration:\n')
cat(sprintf('  base_dir         : %s\n', BASE_DIR))
cat(sprintf('  vcf              : %s\n', VCF_FILE))
cat(sprintf('  env_file         : %s\n', ENV_FILE))
cat(sprintf('  out_rds          : %s\n', OUT_RDS))
cat(sprintf('  n_random_regions : %d\n', n_random_regions))
cat(sprintf('  random_seed      : %d\n', random_seed))
cat(sprintf('  min_panel_snps   : %d\n\n', MIN_PANEL_SNPS_REGION))


# ---- Locus / donor configuration -------------------------------------------

# Chromosome IDs are the underscore form used in VCF after canonicalization.
CHR1 <- 'ENA_OZ004311_OZ004311.1'
CHR2 <- 'ENA_OZ004312_OZ004312.1'
CHR3 <- 'ENA_OZ004313_OZ004313.1'

# Focal regions at the haplotype-block scale, one per locus. These regions
# are BOTH the focal regions for the per-pool test AND the regions excluded
# from the genome-wide panel.
focal_regions <- data.table(
  locus = c('chr1', 'chr2', 'chr3'),
  chrom = c(CHR1,   CHR2,   CHR3),
  start = c(108500000L,  12500000L, 149402438L),
  end   = c(111000000L,  15000000L, 161021387L),
  width = c(  2500000L,   2500000L,  11618950L)   # 2.5 Mb, 2.5 Mb, 11.62 Mb
  # chr3 = CpL3a, DELLY CAL1 breakpoints (paired 3to3/5to5 junction calls; interval bounded by both)
)

# Donor configuration per locus. chr2 carries TWO panel configurations.
locus_donor_config <- list(
  chr1 = list(donor = 'FL1', panel_kind = 'donor_private', env = 'bio05'),
  chr2 = list(donor = 'ME1', panel_kind = 'donor_private', env = 'bio05',
              # second panel: non_FL1_deficit (the FL1>=0.95 fixed, ME1<=0.05, CAL1<=0.05 set;
              # estimator on this panel is 1 - alpha_FL1, which ~= alpha_ME1 + alpha_CAL1)
              secondary = list(donor = 'FL1', panel_kind = 'non_FL1_deficit', env = 'bio05',
                               label = 'chr2_nonFL1_deficit')),
  chr3 = list(donor = 'CAL1', panel_kind = 'donor_private', env = 'pis')
)


# =============================================================================
# Helpers
# =============================================================================

# Given a region's matrices, compute the per-pool alpha (mean AF on the
# subset of region SNPs that are in the donor panel).
#
# `panel_idx` is a logical or integer index into the region's rows; we
# pre-compute panel membership upstream using the reference AF matrix.
#
# Polarization handling: the donor-private criterion (p_D >= 0.95) is defined
# in terms of the ALT allele frequency in the donor reference, computed
# directly from the same VCF in the same allelic orientation. So no flipping
# needed -- panel SNPs are by definition those where ALT is donor-fixed and
# REF-fixed in the others, so the pool ALT frequency directly estimates
# alpha_donor without polarization.
#
# (The allele-frequency cline script polarizes because it works on BayPass
# focal SNPs that aren't pre-filtered for donor-fixity. Here, panel
# construction already imposes the polarization.)
compute_pool_alpha <- function(alt, dp, panel_idx, min_depth = MIN_DEPTH) {
  if (length(panel_idx) == 0L) {
    return(rep(NA_real_, ncol(alt)))
  }
  freq <- alt[panel_idx, , drop = FALSE] / dp[panel_idx, , drop = FALSE]
  freq[dp[panel_idx, , drop = FALSE] < min_depth] <- NA_real_
  colMeans(freq, na.rm = TRUE)
}

# Mean per-pool FST(pool, donor_pool) at region SNPs, simplified pool-level
# estimator: FST_j = (p_i - p_donor)^2 / (mean_p * (1 - mean_p)) at each SNP,
# averaged across SNPs. (This is Reynolds-style per-SNP at the pool level;
# the poolfstat ANOVA estimator is used for the genome-wide pairwise FST.)
compute_mean_fst_to_donor <- function(alt, dp, donor_idx, field_idx,
                                       min_depth = MIN_DEPTH) {
  freq <- alt / dp
  freq[dp < min_depth] <- NA_real_
  p_donor <- freq[, donor_idx]
  
  fst_per_pool <- vapply(field_idx, function(j) {
    p_pool <- freq[, j]
    pbar   <- (p_pool + p_donor) / 2
    num    <- (p_pool - p_donor)^2
    den    <- pbar * (1 - pbar)
    valid  <- is.finite(num) & is.finite(den) & den > 0
    if (sum(valid) == 0L) return(NA_real_)
    mean(num[valid] / den[valid])
  }, numeric(1))
  
  fst_per_pool
}

# Per-region PCA on per-pool AF matrix (field pools only; complete SNPs only).
# Returns: pct_var on PC1 and the PC1 scores vector.
compute_region_pca <- function(alt, dp, field_idx, min_depth = MIN_DEPTH,
                                maf_min = 0.05) {
  freq <- alt[, field_idx] / dp[, field_idx]
  freq[dp[, field_idx] < min_depth] <- NA_real_
  
  # MAF filter (mean over field pools) + complete-rows
  mean_freq <- rowMeans(freq, na.rm = TRUE)
  keep_maf  <- !is.na(mean_freq) & mean_freq >= maf_min & mean_freq <= 1 - maf_min
  freq      <- freq[keep_maf, , drop = FALSE]
  complete  <- rowSums(is.na(freq)) == 0
  freq      <- freq[complete, , drop = FALSE]
  
  if (nrow(freq) < 10L || ncol(freq) < 3L) {
    return(list(pct_var_PC1 = NA_real_, scores_PC1 = rep(NA_real_, length(field_idx)),
                n_snps = nrow(freq)))
  }
  
  centered <- freq - rowMeans(freq)
  pca      <- prcomp(t(centered), center = FALSE, scale. = FALSE)
  pct_var  <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
  
  list(pct_var_PC1 = pct_var[1],
       scores_PC1  = pca$x[, 1],
       n_snps      = nrow(freq))
}

# Test if a candidate random-region interval overlaps any focal broad region.
# A region is rejected (resampled) if it overlaps the union of broad regions
# on any chromosome.
overlaps_broad <- function(chrom, start, end, broad_df) {
  hits <- broad_df$chrom == chrom &
          start <= broad_df$end &
          end   >= broad_df$start
  any(hits)
}


# =============================================================================
# Stage 1: Read VCF (vcfR) and extract pool order + per-pool AD/DP matrices
# =============================================================================
#
# vcfR reads the whole VCF into a single S4 object. From it we extract:
#   - pool order from the @gt column names (skipping the first 'FORMAT' col)
#   - per-pool AD field, parsed as REF,ALT pairs -> alt_all and ref_all
#   - chrom + pos from @fix
#
# Memory: a vcfR object on 11.63M SNPs x 27 pools is ~8-12 GB. Substantially
# heavier than VariantAnnotation's tabix-aware reads, but tolerable on a
# 64 GB cluster node. If RAM becomes a problem we'd switch to bcftools-query
# preprocessing into a flat TSV.

cat('[', format(Sys.time()), '] Stage 1: loading VCF with vcfR (slow step)\n')

vcf <- read.vcfR(VCF_FILE, verbose = FALSE)

# Pool order from the @gt column names (first is FORMAT, rest are samples)
all_pools <- colnames(vcf@gt)[-1]
cat(sprintf('  VCF has %d pools: %s\n',
            length(all_pools), paste(all_pools, collapse = ', ')))

stopifnot(all(REF_POOLS %in% all_pools))
field_pools <- setdiff(all_pools, REF_POOLS)
stopifnot(length(field_pools) == 24L)

field_idx_all <- match(field_pools, all_pools)
ref_idx       <- setNames(match(REF_POOLS, all_pools), REF_POOLS)


# ---- Environmental predictors per field pool --------------------------------

env_dt <- fread(ENV_FILE)
names(env_dt)[names(env_dt) == 'site'] = 'sample'
stopifnot(all(c('sample', 'bio05', 'pis') %in% names(env_dt)))
env_field <- env_dt[sample %in% field_pools]
env_field <- env_field[match(field_pools, sample)]   # reorder to VCF order
stopifnot(identical(env_field$sample, field_pools))


# =============================================================================
# Stage 2: Extract per-pool AD/DP matrices from the vcfR object
# =============================================================================
#
# PoolSNP writes FORMAT=GT:RD:AD:DP:FREQ. We use vcfR::extract.gt to pull
# RD (reference count) and AD (alt count) as numeric matrices; DP = RD + AD.
#
# extract.gt returns a SNPs x samples character matrix; we coerce to numeric.
# Each call costs one full pass through @gt, so we do two calls (one per
# subfield) rather than parsing AD-as-pair (PoolSNP's AD is a scalar, not
# a "REF,ALT" pair like a typical multi-sample VCF).

cat('[', format(Sys.time()), '] Stage 2: extracting AD/RD matrices\n')

alt_all <- extract.gt(vcf, element = 'AD', as.numeric = TRUE)
ref_all <- extract.gt(vcf, element = 'RD', as.numeric = TRUE)
dp_all  <- alt_all + ref_all

stopifnot(identical(colnames(alt_all), all_pools))

snp_chrom <- as.character(vcf@fix[, 'CHROM'])
snp_pos   <- as.integer(vcf@fix[, 'POS'])
n_snps    <- nrow(alt_all)

# Free the parent vcfR object now that we have the matrices + coords
rm(vcf); gc()

cat(sprintf('  Loaded %s SNPs x %d pools\n',
            format(n_snps, big.mark = ','), ncol(alt_all)))

# Reference AFs (with depth filter)
ref_freq <- matrix(NA_real_, nrow = n_snps, ncol = 3,
                   dimnames = list(NULL, REF_POOLS))
for (d in REF_POOLS) {
  i   <- ref_idx[d]
  pf  <- alt_all[, i] / dp_all[, i]
  pf[dp_all[, i] < MIN_DEPTH] <- NA_real_
  ref_freq[, d] <- pf
}

# Panel membership (logical vectors)
panel_membership <- list(
  FL1_private = ref_freq[,'FL1'] >= PANEL_P_DONOR &
                ref_freq[,'ME1'] <= PANEL_P_OTHER &
                ref_freq[,'CAL1'] <= PANEL_P_OTHER,
  ME1_private = ref_freq[,'ME1'] >= PANEL_P_DONOR &
                ref_freq[,'FL1'] <= PANEL_P_OTHER &
                ref_freq[,'CAL1'] <= PANEL_P_OTHER,
  CAL1_private = ref_freq[,'CAL1'] >= PANEL_P_DONOR &
                 ref_freq[,'FL1'] <= PANEL_P_OTHER &
                 ref_freq[,'ME1'] <= PANEL_P_OTHER
)
# non_FL1_deficit: same panel as FL1_private; the alpha estimator differs
# at use time (1 - alpha_FL1 instead of alpha_FL1).
panel_membership$non_FL1_deficit <- panel_membership$FL1_private

# NA-safe: any NA in ref_freq propagates to FALSE
for (p in names(panel_membership)) {
  panel_membership[[p]][is.na(panel_membership[[p]])] <- FALSE
}

cat('  Panel sizes (genome-wide, before broad-region exclusion):\n')
for (p in names(panel_membership)) {
  cat(sprintf('    %-18s : %s SNPs\n', p,
              format(sum(panel_membership[[p]]), big.mark = ',')))
}


# =============================================================================
# Stage 3: Define genome-wide panel index (panel SNPs OUTSIDE focal regions)
# =============================================================================

cat('[', format(Sys.time()), '] Stage 3: defining genome-wide panel\n')

# vcfR returns chromosome names as they appear in the VCF CHROM column,
# which is already underscore form (ENA_OZ004311_OZ004311.1). No
# canonicalization needed -- compare directly.
snp_chrom_underscore <- snp_chrom

in_any_focal <- logical(n_snps)
for (r in seq_len(nrow(focal_regions))) {
  hit <- snp_chrom_underscore == focal_regions$chrom[r] &
         snp_pos >= focal_regions$start[r] &
         snp_pos <= focal_regions$end[r]
  in_any_focal <- in_any_focal | hit
}
cat(sprintf('  SNPs inside any focal region: %s\n',
            format(sum(in_any_focal), big.mark = ',')))
cat(sprintf('  Genome-wide candidate SNPs (focal-excluded): %s\n',
            format(sum(!in_any_focal), big.mark = ',')))


# =============================================================================
# Stage 4: Compute genome-wide alpha per donor (the alpha_genome baseline)
# =============================================================================

cat('[', format(Sys.time()), '] Stage 4: computing genome-wide alpha per donor\n')

# Genome-wide alpha = mean over panel SNPs (broad-excluded) of pool ALT
# frequency. We compute once per donor configuration; the result is a
# 24-length vector keyed by field pool name.

compute_genome_alpha <- function(panel_name, sign = 1) {
  # sign = +1 for donor_private (alpha = mean ALT)
  # sign = -1 for non_FL1_deficit (alpha = 1 - mean ALT_FL1_panel)
  panel_idx <- which(panel_membership[[panel_name]] & !in_any_focal)
  if (length(panel_idx) == 0L) stop('Empty genome-wide panel for ', panel_name)
  
  freq_field <- alt_all[panel_idx, field_idx_all] /
                dp_all [panel_idx, field_idx_all]
  freq_field[dp_all[panel_idx, field_idx_all] < MIN_DEPTH] <- NA_real_
  
  alpha <- colMeans(freq_field, na.rm = TRUE)
  if (sign == -1) alpha <- 1 - alpha
  names(alpha) <- field_pools
  alpha
}

genome_panel_alpha <- list(
  FL1_private      = compute_genome_alpha('FL1_private',  sign =  1),
  ME1_private      = compute_genome_alpha('ME1_private',  sign =  1),
  CAL1_private     = compute_genome_alpha('CAL1_private', sign =  1),
  non_FL1_deficit  = compute_genome_alpha('non_FL1_deficit', sign = -1)
)

cat('  Genome-wide alpha summary (per donor; mean across 24 field pools):\n')
for (p in names(genome_panel_alpha)) {
  a <- genome_panel_alpha[[p]]
  cat(sprintf('    %-18s : mean=%.3f  sd=%.3f  range=[%.3f, %.3f]\n',
              p, mean(a, na.rm = TRUE), sd(a, na.rm = TRUE),
              min(a, na.rm = TRUE), max(a, na.rm = TRUE)))
}


# =============================================================================
# Stage 5: Build null distributions per locus by sampling random regions
# =============================================================================

cat('[', format(Sys.time()), '] Stage 5: sampling random regions\n')

# Chromosome sizes for region sampling, from the idCulPipi1 assembly (Hesson
# et al. 2025): chr1 = 128.5 Mb, chr2 = 213.1 Mb, chr3 = 190.9 Mb.
# (We use only autosomes; the three OZ contigs ARE the three autosomes.)
chrom_sizes <- c(
  ENA_OZ004311_OZ004311.1 = 128500000L,
  ENA_OZ004312_OZ004312.1 = 213100000L,
  ENA_OZ004313_OZ004313.1 = 190900000L
)
genome_total <- sum(chrom_sizes)

# Sample N random regions of a given width across autosomes, weighting by
# chromosome size. Rejection-sample to exclude broad regions.
sample_random_regions <- function(n, width, broad_df, seed) {
  set.seed(seed)
  out <- vector('list', n)
  attempts <- 0L
  max_attempts <- n * 20L   # safety: avoid infinite loop
  i <- 1L
  while (i <= n && attempts < max_attempts) {
    attempts <- attempts + 1L
    chrom <- sample(names(chrom_sizes), 1L,
                    prob = chrom_sizes / genome_total)
    max_start <- chrom_sizes[chrom] - width
    if (max_start <= 1L) next
    start <- sample.int(max_start, 1L)
    end   <- start + width - 1L
    if (overlaps_broad(chrom, start, end, broad_df)) next
    out[[i]] <- list(chrom = chrom, start = start, end = end)
    i <- i + 1L
  }
  if (i <= n) {
    warning(sprintf('Only sampled %d of %d regions before hitting attempt cap', i - 1L, n))
    out <- out[seq_len(i - 1L)]
  }
  rbindlist(lapply(out, as.data.table))
}

# For each locus, compute null statistics on each random region
compute_null_for_locus <- function(locus_name, donor, panel_name, env_var,
                                    width, alpha_sign, regions_dt, 
                                    label = locus_name) {
  cat(sprintf('  [%s -- panel %s -- width %d bp -- %d regions]\n',
              label, panel_name, width, nrow(regions_dt)))
  
  alpha_genome <- genome_panel_alpha[[panel_name]]
  env_vec      <- env_field[[env_var]]
  donor_pool_idx <- ref_idx[donor]
  
  results <- vector('list', nrow(regions_dt))
  start_time <- Sys.time()
  
  for (r in seq_len(nrow(regions_dt))) {
    rr <- regions_dt[r]
    
    # Pull SNPs in region (chrom in underscore form for matching)
    in_region <- snp_chrom_underscore == rr$chrom &
                 snp_pos >= rr$start &
                 snp_pos <= rr$end
    if (sum(in_region) < MIN_PANEL_SNPS_REGION) {
      results[[r]] <- data.table(
        chrom = rr$chrom, start = rr$start, end = rr$end,
        n_snps_in_region = sum(in_region), n_panel_snps = 0L,
        delta_marginal = NA_real_, delta_slope = NA_real_,
        pc1_pct_var = NA_real_, pc1_fst_cor = NA_real_,
        fst_donor_var = NA_real_
      )
      next
    }
    
    region_idx <- which(in_region)
    
    # ---- Delta statistics (panel-restricted) -----------------------------
    
    panel_in_region <- which(panel_membership[[panel_name]][region_idx])
    panel_region_idx <- region_idx[panel_in_region]
    
    if (length(panel_region_idx) >= MIN_PANEL_SNPS_REGION) {
      freq_field <- alt_all[panel_region_idx, field_idx_all, drop = FALSE] /
                    dp_all [panel_region_idx, field_idx_all, drop = FALSE]
      freq_field[dp_all[panel_region_idx, field_idx_all, drop = FALSE] < MIN_DEPTH] <- NA_real_
      alpha_region <- colMeans(freq_field, na.rm = TRUE)
      if (alpha_sign == -1) alpha_region <- 1 - alpha_region
      
      delta_i        <- alpha_region - alpha_genome
      delta_marginal <- median(delta_i, na.rm = TRUE)
      slope_fit      <- tryCatch(coef(lm(delta_i ~ env_vec))[2], error = function(e) NA_real_)
      delta_slope    <- as.numeric(slope_fit)
    } else {
      delta_marginal <- NA_real_
      delta_slope    <- NA_real_
      panel_region_idx <- integer(0)
    }
    
    # ---- PCA-FST statistics (all region SNPs, not panel-restricted) ------
    
    pca_res <- compute_region_pca(
      alt_all[region_idx, , drop = FALSE],
      dp_all [region_idx, , drop = FALSE],
      field_idx_all
    )
    
    fst_per_pool <- compute_mean_fst_to_donor(
      alt_all[region_idx, , drop = FALSE],
      dp_all [region_idx, , drop = FALSE],
      donor_idx = donor_pool_idx,
      field_idx = field_idx_all
    )
    
    pc1_fst_cor   <- if (all(is.finite(pca_res$scores_PC1)) && all(is.finite(fst_per_pool)))
                       abs(cor(pca_res$scores_PC1, fst_per_pool)) else NA_real_
    fst_donor_var <- var(fst_per_pool, na.rm = TRUE)
    
    results[[r]] <- data.table(
      chrom = rr$chrom, start = rr$start, end = rr$end,
      n_snps_in_region = length(region_idx),
      n_panel_snps     = length(panel_region_idx),
      delta_marginal   = delta_marginal,
      delta_slope      = delta_slope,
      pc1_pct_var      = pca_res$pct_var_PC1,
      pc1_fst_cor      = pc1_fst_cor,
      fst_donor_var    = fst_donor_var
    )
    
    if (r %% 50L == 0L || r == nrow(regions_dt)) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = 'mins'))
      cat(sprintf('    [%d / %d]  elapsed=%.1f min\n', r, nrow(regions_dt), elapsed))
    }
  }
  
  rbindlist(results)
}


# ---- Per-locus null builds --------------------------------------------------

null_list <- list()

for (locus_name in c('chr1', 'chr2', 'chr3')) {
  cfg <- locus_donor_config[[locus_name]]
  width <- focal_regions$width[focal_regions$locus == locus_name]
  
  # Sample random regions once for this locus width; reuse across panels
  cat(sprintf('\n--- locus %s (donor=%s, width=%d bp) ---\n',
              locus_name, cfg$donor, width))
  regions_dt <- sample_random_regions(
    n        = n_random_regions,
    width    = width,
    broad_df = focal_regions,
    seed     = random_seed + match(locus_name, c('chr1','chr2','chr3'))
  )
  
  # Primary panel
  primary_panel <- paste0(cfg$donor, '_private')
  null_primary  <- compute_null_for_locus(
    locus_name  = locus_name,
    donor       = cfg$donor,
    panel_name  = primary_panel,
    env_var     = cfg$env,
    width       = width,
    alpha_sign  = 1,
    regions_dt  = regions_dt,
    label       = paste0(locus_name, '_', primary_panel)
  )
  null_list[[locus_name]] <- list(primary = null_primary)
  
  # Secondary panel (chr2 only)
  if (!is.null(cfg$secondary)) {
    sec <- cfg$secondary
    null_secondary <- compute_null_for_locus(
      locus_name = locus_name,
      donor      = sec$donor,
      panel_name = sec$panel_kind,
      env_var    = sec$env,
      width      = width,
      alpha_sign = -1,                           # 1 - alpha_FL1
      regions_dt = regions_dt,
      label      = sec$label
    )
    null_list[[locus_name]]$secondary <- null_secondary
  }
}


# =============================================================================
# Stage 6: Save artifact
# =============================================================================

out <- list(
  config = list(
    n_random_regions = n_random_regions,
    random_seed      = random_seed,
    panel_p_donor    = PANEL_P_DONOR,
    panel_p_other    = PANEL_P_OTHER,
    min_depth        = MIN_DEPTH,
    min_panel_snps_region = MIN_PANEL_SNPS_REGION,
    locus_regions    = focal_regions,
    field_pools      = field_pools,
    timestamp        = format(Sys.time())
  ),
  panels             = lapply(panel_membership, sum),     # just counts; raw bool vecs are big
  panel_membership   = panel_membership,                  # full vectors for downstream reuse
  genome_panel_alpha = genome_panel_alpha,
  null               = null_list
)

saveRDS(out, OUT_RDS)
cat(sprintf('\n[%s] Wrote: %s\n', format(Sys.time()), OUT_RDS))
cat(sprintf('         File size: %.1f MB\n', file.size(OUT_RDS) / 1e6))

# ---- Summary diagnostics ----------------------------------------------------

cat('\n=== Null distribution summaries (median +/- IQR) ===\n')
for (l in names(null_list)) {
  for (k in names(null_list[[l]])) {
    n <- null_list[[l]][[k]]
    cat(sprintf('  %s / %s  (n_regions valid = %d)\n', l, k, sum(!is.na(n$delta_marginal))))
    for (stat in c('delta_marginal','delta_slope','pc1_pct_var','pc1_fst_cor','fst_donor_var')) {
      v <- n[[stat]]
      if (all(is.na(v))) next
      cat(sprintf('    %-16s median=%.4f  IQR=[%.4f, %.4f]\n', stat,
                  median(v, na.rm = TRUE),
                  quantile(v, 0.25, na.rm = TRUE),
                  quantile(v, 0.75, na.rm = TRUE)))
    }
  }
}

cat('\nDone.\n')
