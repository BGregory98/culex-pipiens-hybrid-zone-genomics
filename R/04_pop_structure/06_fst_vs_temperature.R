# =============================================================================
# R/04_pop_structure/06_fst_vs_temperature.R
# =============================================================================
# Genome-wide FST between each field pool and each reference pool against
# maximum summer temperature (bio05), before (A) and after (B) controlling for
# genome-wide expected heterozygosity (He).
#
# (A) FST ~ bio05 per reference, OLS fit with 95% band.
# (B) Added-variable plots from FST ~ bio05 + He: residuals of FST ~ He against
#     residuals of bio05 ~ He. The fitted slope equals the bio05 coefficient of
#     the joint model, and its P value is the one reported.
#
# Manuscript: Figure S5; Results (P = 0.72 / 0.006 / < 0.001); Supplementary
#             Results "Genome-wide diversity and differentiation from the
#             reference pools"
#
# Inputs:  data/genome_wide_FST_for_field_vs_references.csv  (Table S2; R/04_pop_structure/01_pairwise_fst.R)
#          data/site_coords_with_spatial_data.csv            (bio05; R/00_setup/01_site_predictors.R)
#          results/multilocus/Table_per_pool_diversity.csv   (He; R/04_pop_structure/05_pool_diversity.R)
# Outputs: figures/final/supplement/Figure_S5.png
#          results/multilocus/Table_fst_vs_bio05_partial.csv (per-reference statistics)
#
# Runtime: < 1 min
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/ref_colors.R")   # ref_labels

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

FST_FILE <- "data/genome_wide_FST_for_field_vs_references.csv"
ENV_FILE <- "data/site_coords_with_spatial_data.csv"
HET_FILE <- "results/multilocus/Table_per_pool_diversity.csv"
need_files(FST_FILE, ENV_FILE, HET_FILE)

FIG_OUT   <- "figures/final/supplement/Figure_S5.png"
TABLE_OUT <- "results/multilocus/Table_fst_vs_bio05_partial.csv"
dir.create(dirname(FIG_OUT), recursive = TRUE, showWarnings = FALSE)

fst <- fread(FST_FILE)
env <- fread(ENV_FILE)[, .(sample = site, bio05)]
het <- fread(HET_FILE)[, .(sample = pool, He = exp_het)]

d <- merge(merge(fst, env, by = "sample"), het, by = "sample")
d[, ref := factor(ref_labels[ref], levels = c("PIP", "MOL", "QUI"))]
stopifnot(d[, .N, by = ref][, all(N == 24)])

# ---- statistics --------------------------------------------------------------
stats <- d[, {
  m1 <- summary(lm(gw_fst ~ bio05))
  m2 <- summary(lm(gw_fst ~ bio05 + He))
  .(slope_raw = coef(m1)["bio05", 1], adj_r2_raw = m1$adj.r.squared, p_raw = coef(m1)["bio05", 4],
    slope_partial = coef(m2)["bio05", 1], p_partial = coef(m2)["bio05", 4])
}, by = ref]
print(stats, digits = 3)
fwrite(stats, TABLE_OUT)

# ---- added-variable coordinates ---------------------------------------------
d[, `:=`(fst_res = resid(lm(gw_fst ~ He)), bio05_res = resid(lm(bio05 ~ He))), by = ref]

theme_s5 <- theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey85"))

pA <- ggplot(d, aes(bio05, gw_fst)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#3366FF") +
  geom_point(size = 1.8) +
  facet_wrap(~ ref, ncol = 1, scales = "free_y") +
  labs(x = "Maximum summer temperature (°C)", y = expression(F[ST])) +
  theme_s5

pB <- ggplot(d, aes(bio05_res, fst_res)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#3366FF") +
  geom_point(size = 1.8) +
  facet_wrap(~ ref, ncol = 1, scales = "free_y") +
  labs(x = "Temperature residual (°C)", y = expression(F[ST]~"residual")) +
  theme_s5

fig <- pA + pB + plot_annotation(tag_levels = "A")
ggsave(FIG_OUT, fig, width = 8, height = 9, dpi = 300, bg = "white")
cat("wrote", FIG_OUT, "\n")
