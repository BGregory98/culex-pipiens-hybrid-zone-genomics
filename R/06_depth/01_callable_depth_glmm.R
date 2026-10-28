# =============================================================================
# R/06_depth/01_callable_depth_glmm.R
# =============================================================================
# Callable-depth GLMM for the CpL1a CYP6 paralogs and the CpL2a Yan ortholog.
# Per-gene total CDS depth at PoolSNP-callable positions (13,741 genes x 24
# field pools) is modelled in glmmTMB with log genome-wide callable depth as an
# offset and a per-gene random intercept. M0 fits z-scored bio05 genome-wide
# (Poisson, refit as nbinom2 when overdispersed); M1, fit once per locus, adds
# a focal-gene indicator and its interaction with bio05, so the interaction
# tests the focal genes against the genome-wide depth-temperature slope.
# Draws Figure 6: per-position depth ratio at each CYP6 paralog against bio05.
#
# Manuscript: Figure 6; Results "Cx. quinquefasciatus ancestry and sequencing
#             read depth at a CYP6 cluster rise with summer temperature"
#             (18 times the genome-wide slope, ~31% more depth) and CpL2a
#             (Yan, P = 0.60); Methods and Supplementary Methods
#             "Callable-position read-depth analysis" (Poisson dispersion)
#
# Inputs:  data/glmm_v4/cds_depth_per_gene_pool.tsv         (cluster/06_depth/glmm/04_aggregate_cds_depth.R)
#          data/glmm_v4/genome_callable_depth_per_pool.tsv  (cluster/06_depth/glmm/04_aggregate_cds_depth.R)
#          data/site_coords_with_spatial_data.csv           (R/00_setup/01_site_predictors.R)
#          data/sample_categories.csv                       (sample metadata)
# Outputs: results/glmm/glmm_model_summary.txt   model summaries (M0, M1 per locus)
#          results/glmm/glmm_focal_effects.tsv   focal-gene interaction estimate per locus
#          figures/final/Figure_6.png
#
# Runtime: a few minutes (three GLMM fits on ~330,000 gene x pool rows)
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/sample_metadata.R")   # field_names
source("R/helpers/loci.R")              # CYP6_GENES, YAN_GENE

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
})


# =============================================================================
# Focal-gene set per locus
# =============================================================================
# chr1: the six CpL1a CYP6 paralogs; chr2: the CpL2a Yan ortholog.

focal_genes <- list(
  chr1 = CYP6_GENES,
  chr2 = YAN_GENE
)

# P450Atlas assignments of the six paralogs, for the Figure 6 legend.
CYP6_LABELS <- c(CQUJHB009005 = 'CYP6AA7-like', CQUJHB008986 = 'CYP6AA8-like',
                 CQUJHB009148 = 'CYP6AA9-like', CQUJHB017007 = 'CYP6BZ2-like',
                 CQUJHB001163 = 'CYP6CC2-like', CQUJHB007871 = 'CYP6BB4-like')
stopifnot(setequal(names(CYP6_LABELS), CYP6_GENES))

# =============================================================================
# INPUTS
# =============================================================================

GLMM_DATA_DIR <- 'data/glmm_v4'
GLMM_OUT_DIR  <- 'results/glmm'
FIG_OUT       <- 'figures/final/Figure_6.png'

CDS_FILE  <- file.path(GLMM_DATA_DIR, 'cds_depth_per_gene_pool.tsv')
GW_FILE   <- file.path(GLMM_DATA_DIR, 'genome_callable_depth_per_pool.tsv')
SITE_FILE <- 'data/site_coords_with_spatial_data.csv'
need_files(CDS_FILE, GW_FILE, SITE_FILE)

dir.create(GLMM_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(FIG_OUT), recursive = TRUE, showWarnings = FALSE)

cds_depth <- fread(CDS_FILE)
gw_depth  <- fread(GW_FILE)

stopifnot(c('gene_id', 'pool', 'total_depth', 'n_callable_pos') %in% names(cds_depth))
stopifnot(c('pool', 'genome_callable_depth') %in% names(gw_depth))

cat(sprintf('[INFO] cds_depth: %s rows, %s unique gene_ids, %s pools\n',
            format(nrow(cds_depth), big.mark = ','),
            format(uniqueN(cds_depth$gene_id), big.mark = ','),
            uniqueN(cds_depth$pool)))


# =============================================================================
# Environmental data: bio05 per field pool
# =============================================================================
# bio05 as the mean within a 500 m buffer around each site (Table S1).

env_data <- read.csv(SITE_FILE) %>%
  dplyr::filter(site %in% field_names) %>%
  dplyr::rename(pool = site) %>%
  dplyr::select(pool, bio05)


# =============================================================================
# Assemble the modeling table
# =============================================================================

model_data <- cds_depth %>%
  left_join(gw_depth, by = 'pool') %>%
  left_join(env_data, by = 'pool') %>%
  mutate(
    # Locus assignment: chr1 if gene is in chr1 focal, chr2 if in chr2 focal,
    # NA if neither. The genome-wide null is the entire set of genes; we use
    # is_focal indicators per locus rather than dropping non-focal genes.
    is_focal_chr1 = gene_id %in% focal_genes$chr1,
    is_focal_chr2 = gene_id %in% focal_genes$chr2,
    log_gw        = log(genome_callable_depth),
    bio05_z       = as.numeric(scale(bio05))   # z-score for numerical stability
  )

# Sanity check: every focal gene present in the data?
for (locus in names(focal_genes)) {
  expected <- focal_genes[[locus]]
  present  <- intersect(expected, model_data$gene_id)
  missing  <- setdiff(expected, model_data$gene_id)
  cat(sprintf('[INFO] %s: %d / %d focal genes present in model data\n',
              locus, length(present), length(expected)))
  if (length(missing) > 0) {
    cat(sprintf('         missing: %s\n', paste(missing, collapse = ', ')))
    cat(sprintf('         (typically due to the minimum-positions filter in the aggregation step dropping low-coverage genes)\n'))
  }
}


# =============================================================================
# Fit M0 -- baseline, genome-wide
# =============================================================================
# total_depth ~ bio05 + offset(log(genome_callable_depth)) + (1 | gene_id)
#
# Start with Poisson, refit with nbinom2 if overdispersed.

cat('\n========== M0 (baseline, genome-wide) ==========\n')

fit_m0_pois <- glmmTMB(
  total_depth ~ bio05_z + offset(log_gw) + (1 | gene_id),
  data   = model_data,
  family = poisson()
)

# Check Poisson overdispersion. Pearson residual variance / df_resid should be ~ 1.
pearson_resid <- residuals(fit_m0_pois, type = 'pearson')
df_resid      <- df.residual(fit_m0_pois)
dispersion    <- sum(pearson_resid^2) / df_resid
cat(sprintf('[INFO] M0 Poisson dispersion: %.2f (target ~ 1)\n', dispersion))

if (dispersion > 2) {
  cat('[INFO] Poisson overdispersed; refitting with nbinom2\n')
  fit_m0 <- glmmTMB(
    total_depth ~ bio05_z + offset(log_gw) + (1 | gene_id),
    data   = model_data,
    family = nbinom2()
  )
  m0_family <- 'nbinom2'
} else {
  fit_m0    <- fit_m0_pois
  m0_family <- 'poisson'
}

m0_summary <- summary(fit_m0)
cat(sprintf('M0 family: %s\n', m0_family))
print(m0_summary$coefficients$cond)


# =============================================================================
# Fit M1 -- focal interaction, per locus
# =============================================================================
# Per-locus M1: fit separately for chr1 and chr2 rather than one combined model.
# Reason: a combined model averages signals across loci, which could hide a
# strong signal at one locus and a null at the other.

cat('\n========== M1 (focal interaction, per locus) ==========\n')

fit_m1 <- list()
m1_focal_effects <- list()

for (locus in names(focal_genes)) {
  cat(sprintf('\n---- %s ----\n', locus))

  is_focal_col <- paste0('is_focal_', locus)
  model_data$is_focal_current <- model_data[[is_focal_col]]

  fit_m1_locus <- glmmTMB(
    total_depth ~ bio05_z * is_focal_current + offset(log_gw) + (1 | gene_id),
    data   = model_data,
    family = if (m0_family == 'nbinom2') nbinom2() else poisson()
  )

  fit_m1[[locus]] <- fit_m1_locus
  s <- summary(fit_m1_locus)
  print(s$coefficients$cond)

  # Extract the bio05_z:is_focal_currentTRUE coefficient (the focal-interaction effect).
  coefs <- s$coefficients$cond
  interaction_row <- grep('bio05_z:is_focal', rownames(coefs))
  if (length(interaction_row) == 1) {
    m1_focal_effects[[locus]] <- data.frame(
      locus     = locus,
      n_focal   = sum(model_data$is_focal_current),
      n_focal_genes = length(focal_genes[[locus]]),
      estimate  = coefs[interaction_row, 'Estimate'],
      std_error = coefs[interaction_row, 'Std. Error'],
      z_value   = coefs[interaction_row, 'z value'],
      p_value   = coefs[interaction_row, 'Pr(>|z|)']
    )
  }
}

m1_effects_df <- do.call(rbind, m1_focal_effects)


# =============================================================================
# Write text summaries
# =============================================================================

summary_path <- file.path(GLMM_OUT_DIR, 'glmm_model_summary.txt')
sink(summary_path)
cat('Callable-depth GLMM summaries\n')
cat(sprintf('Generated: %s\n', Sys.time()))
cat(sprintf('Model family: %s\n\n', m0_family))

cat('==================== M0 (genome-wide) ====================\n')
print(summary(fit_m0))

for (locus in names(fit_m1)) {
  cat(sprintf('\n==================== M1 (%s focal interaction) ====================\n', locus))
  print(summary(fit_m1[[locus]]))
}

cat('\n==================== Focal interaction effects (M1 per locus) ====================\n')
print(m1_effects_df, row.names = FALSE)
sink()

fwrite(m1_effects_df, file.path(GLMM_OUT_DIR, 'glmm_focal_effects.tsv'), sep = '\t')
cat(sprintf('\n[INFO] wrote model summaries to %s\n', summary_path))


# =============================================================================
# Figure 6
# =============================================================================
# For each CYP6 paralog, per-pool per-position depth at the gene, normalized
# by the pool's genome-wide per-position average (y = 1.0: gene depth matches
# genome-average depth; y = 1.5: 1.5x the genome-average per-position depth).
#
# Note: the GLMM models total_depth (the per-gene-summed depth) with an
# offset for genome_callable_depth; this visualization shows the same
# information rescaled to per-position units for readability. The per-gene
# random intercept in the GLMM absorbs gene-length variation, so the model
# coefficients describe the same relationship as the per-position ratio
# plotted here.

# Total callable positions genome-wide, from the callable-depth extraction on
# the cluster (positions passing min-cov / max-cov). Used as the denominator
# for the genome-average per-position depth.
GLOBAL_CALLABLE_POSITIONS <- 155123545L

make_locus_panel <- function(locus, model_data, focal_genes) {

  focal_ids <- focal_genes[[locus]]

  plot_data <- model_data %>%
    mutate(
      # Per-base depth at the gene in this pool
      gene_per_pos_depth   = total_depth / n_callable_pos,
      # Per-base depth genome-wide in this pool
      gw_per_pos_depth     = genome_callable_depth / GLOBAL_CALLABLE_POSITIONS,
      # Ratio: 1.0 = gene matches genome average; >1 = elevated
      depth_ratio          = gene_per_pos_depth / gw_per_pos_depth,
      gene_type            = ifelse(gene_id %in% focal_ids, 'focal', 'genome-wide background')
    )

  # Focal-gene depth ratios per pool, plus background distribution per pool for context.
  focal_data <- plot_data %>% filter(gene_id %in% focal_ids)
  bg_summary <- plot_data %>%
    filter(!(gene_id %in% focal_ids)) %>%
    group_by(pool, bio05) %>%
    summarize(bg_median = median(depth_ratio, na.rm = TRUE),
              bg_q25    = quantile(depth_ratio, 0.25, na.rm = TRUE),
              bg_q75    = quantile(depth_ratio, 0.75, na.rm = TRUE),
              .groups   = 'drop')

  ggplot() +
    # Reference line at y = 1 (genome-average per-position depth)
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey60', linewidth = 0.4) +
    # Background IQR ribbon
    geom_ribbon(data = bg_summary,
                aes(x = bio05, ymin = bg_q25, ymax = bg_q75),
                fill = 'grey80', alpha = 0.5) +
    geom_line(data = bg_summary,
              aes(x = bio05, y = bg_median),
              color = 'grey50', linewidth = 0.8) +
    # Focal genes: one line per gene, points per pool
    geom_point(data = focal_data,
               aes(x = bio05, y = depth_ratio, color = gene_id),
               size = 1.8, alpha = 0.4) +
    geom_smooth(data = focal_data,
                aes(x = bio05, y = depth_ratio, color = gene_id, group = gene_id),
                method = 'lm', se = FALSE, linewidth = 0.7) +
    scale_color_brewer(palette = 'Dark2', name = NULL,
                       labels = function(ids) unname(CYP6_LABELS[ids])) +
    labs(
      x     = 'Maximum Summer Temperature (°C)',
      y     = 'Per-position depth at gene\n/ genome-wide per-position depth'    ) +
    theme_bw(base_size = 18) +
    theme(legend.position = 'right')
}

fig6 <- make_locus_panel('chr1', model_data, focal_genes)

ggsave(FIG_OUT, fig6, width = 12, height = 9,
       dpi = 300, bg = 'white')

cat(sprintf('[INFO] wrote figure to %s\n', FIG_OUT))


# =============================================================================
# Console summary
# =============================================================================

cat('\nM0 bio05 fixed-effect (genome-wide):\n')
m0_coefs <- summary(fit_m0)$coefficients$cond
print(m0_coefs[grep('bio05', rownames(m0_coefs)), , drop = FALSE])

cat('\nM1 bio05:is_focal interaction (focal genes vs genome-wide background):\n')
print(m1_effects_df, row.names = FALSE)

beta_int  <- fixef(fit_m1$chr1)$cond["bio05_z:is_focal_currentTRUE"]
beta_bg   <- fixef(fit_m1$chr1)$cond["bio05_z"]
bio05_z_span <- diff(range(model_data$bio05_z))
bio05_span_C <- diff(range(model_data$bio05))   # for °C interpretation

cat(sprintf("CpL1a focal interaction / genome-wide bio05 slope: %.1f\n",
            beta_int / beta_bg))
cat(sprintf("Empirical bio05_z range: %.2f SDs (%.2f °C)\n",
            bio05_z_span, bio05_span_C))
cat(sprintf("Focal excess across this range: %+.1f%%\n",
            100 * (exp(bio05_z_span * beta_int) - 1)))
cat(sprintf("Per-°C focal excess:           %+.1f%% per °C\n",
            100 * (exp(beta_int / (bio05_span_C / bio05_z_span)) - 1)))

cat('\n[OK] done.\n')
