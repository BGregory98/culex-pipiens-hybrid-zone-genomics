# =============================================================================
# R/07_multi_locus/17_spatial_covariate_test.R
# =============================================================================
# Whether each private-panel ancestry-environment association survives adding
# latitude and longitude as covariates. No shuffle of bio05 across the 24 sites
# reproduces its spatial autocorrelation (16_spatial_permutation.R), so
# permutation inference on it is anti-conservative; the spatial structure is
# modelled as a covariate instead. Per-pool panel alpha is regressed on the
# locus predictor with and without lat + lon, and each predictor is regressed
# on lat + lon to show how much of it the coordinates explain.
#
# Manuscript: Supplementary Results "Temperature is collinear with position along
#             the corridor" (coordinate control on raw panel ancestry);
#             Supplementary Methods "Geographic and metropolitan-area controls"
#
# Inputs:  results/multilocus/Table_private_panel_alpha_per_pool.csv  (R/07_multi_locus/12_private_panel_clines.R)
#          data/site_coords_with_spatial_data_allvars.csv             (R/00_setup/01_site_predictors.R)
# Outputs: none (console only)
#
# Runtime: seconds
# =============================================================================

source("R/helpers/config.R")
suppressMessages(library(data.table))

need_files("results/multilocus/Table_private_panel_alpha_per_pool.csv",
           "data/site_coords_with_spatial_data_allvars.csv")

alpha <- fread("results/multilocus/Table_private_panel_alpha_per_pool.csv", showProgress = FALSE)
env   <- fread("data/site_coords_with_spatial_data_allvars.csv", showProgress = FALSE)
sites <- sort(unique(alpha$pool))
setkey(env, site); e <- env[sites]

PAIRS <- list(
  list(loc = "CpL1a", pn = "QUI-private", v = "bio05"),
  list(loc = "CpL1a", pn = "MOL-private", v = "bio05"),
  list(loc = "CpL2a", pn = "QUI-private", v = "bio05"),
  list(loc = "CpL2a", pn = "MOL-private", v = "bio05"),
  list(loc = "CpL3a", pn = "MOL-private", v = "pis")
)

say("%-6s %-12s %-6s %11s %11s %11s %11s", "locus", "panel", "env",
    "P(env only)", "P(+lat/lon)", "adjR2 only", "adjR2 +geo")
for (P in PAIRS) {
  d <- alpha[locus == P$loc & panel == P$pn]
  setkey(d, pool); d <- d[sites]
  if (!nrow(d) || sd(d$alpha) == 0) next
  x <- as.numeric(e[[P$v]])
  f1 <- lm(d$alpha ~ x)
  f2 <- lm(d$alpha ~ x + e$lat + e$lon)
  say("%-6s %-12s %-6s %11.4g %11.4g %11.3f %11.3f",
      P$loc, P$pn, P$v,
      summary(f1)$coefficients["x", 4],
      summary(f2)$coefficients["x", 4],
      summary(f1)$adj.r.squared, summary(f2)$adj.r.squared)
}

say("")
say("=== for reference: how much of each predictor geography explains ===")
for (v in c("bio05", "pis")) {
  s <- summary(lm(as.numeric(e[[v]]) ~ e$lat + e$lon))
  say("  %-6s ~ lat + lon : R2 = %.3f", v, s$r.squared)
}
say("")
say("Reading: where P(+lat/lon) stays small, the environmental association is not")
say("reducible to position along the corridor. Where it inflates, it is.")
