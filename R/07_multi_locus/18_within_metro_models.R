# =============================================================================
# R/07_multi_locus/18_within_metro_models.R
# =============================================================================
# Separates the temperature gradient from the Baltimore-Washington contrast.
# Temperature varies chiefly between the cooler Baltimore and the warmer
# Washington cores, while impervious surface varies within each, so this asks
# what survives WITHIN metropolitan area. Sites are assigned to the nearer city
# centre by great-circle distance.
#   1. For each predictor: variance explained by coordinates and by metro, its
#      total and within-metro SD, and its correlation with bio05.
#   2. For each locus Delta (per-population ancestry test): P for each predictor
#      alone, with a metro term, and with latitude + longitude; bio05 fitted
#      within each metro separately; and bio05 + impervious surface within
#      Washington. Section 1 also gives the within-metro bio05-impervious
#      correlations and bio05's R2 on coordinates within each metro.
#   3. Delta ~ bio05 + impervious surface + metro, for every locus and panel.
#
# Manuscript: Supplementary Results "Temperature is collinear with position along
#             the corridor" (R2 on coordinates; within-area fits and bio05 spans;
#             coordinate-controlled P on Delta; both axes with a metro term);
#             Supplementary Methods "Geographic and metropolitan-area controls"
#             (11 and 13 sites); Results CpL1a and CpL3a sections (metropolitan-
#             area term)
#
# Inputs:  results/multilocus/per_pool_ancestry_per_pool.csv  (R/07_multi_locus/01_per_pool_ancestry.R)
#          data/site_coords_with_spatial_data_allvars.csv     (R/00_setup/01_site_predictors.R)
#          data/sample_categories.csv
# Outputs: none (console only)
#
# Runtime: seconds
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/sample_metadata.R")
source("R/helpers/geo.R")
suppressMessages(library(data.table))

need_files("results/multilocus/per_pool_ancestry_per_pool.csv",
           "data/site_coords_with_spatial_data_allvars.csv")

d    <- fread("results/multilocus/per_pool_ancestry_per_pool.csv", showProgress = FALSE)
env  <- fread("data/site_coords_with_spatial_data_allvars.csv", showProgress = FALSE)
setkey(env, site); e <- env[field_names]
e$class <- sample_info_all$sample_type_2[match(e$site, sample_info_all$sample)]

# ---- metro assignment: nearer of Baltimore / Washington city centres ----
e$d_bal <- hav(e$lat, e$lon, BAL[1], BAL[2])
e$d_dc  <- hav(e$lat, e$lon, DC[1],  DC[2])
e$metro <- ifelse(e$d_bal < e$d_dc, "Baltimore", "Washington")
say("metro split: Baltimore %d, Washington %d",
    sum(e$metro == "Baltimore"), sum(e$metro == "Washington"))
print(e[, .(site, metro, class, bio05 = round(bio05, 2), pis = round(pis, 1))],
      row.names = FALSE)

# ---- 1. variance partition of each predictor ----
COVS <- c("bio05", "bio07", "bio18", "elevation", "pis", "tree_canopy", "pd", "wtd")
say("")
say("=== how confounded is each predictor? ===")
say("%-12s %9s %9s %10s %10s %9s", "predictor", "R2 lat/lon", "R2 metro",
    "SD total", "SD within", "r w/bio05")
for (v in COVS) {
  x <- as.numeric(e[[v]])
  r_geo   <- summary(lm(x ~ e$lat + e$lon))$r.squared
  r_metro <- summary(lm(x ~ e$metro))$r.squared
  sd_w    <- sqrt(mean(tapply(x, e$metro, function(z) var(z) * (length(z) - 1) /
                              length(z))))
  say("%-12s %9.3f %9.3f %10.4g %10.4g %9.3f",
      v, r_geo, r_metro, sd(x), sd_w, cor(x, as.numeric(e$bio05)))
}
say("")
say("  bio05 within-metro range: Baltimore %.2f-%.2f C, Washington %.2f-%.2f C",
    min(e$bio05[e$metro == "Baltimore"]), max(e$bio05[e$metro == "Baltimore"]),
    min(e$bio05[e$metro == "Washington"]), max(e$bio05[e$metro == "Washington"]))
for (mm in c("Washington", "Baltimore")) {
  k  <- e$metro == mm
  ct <- cor.test(as.numeric(e$bio05[k]), as.numeric(e$pis[k]))
  say("  %-10s bio05 vs impervious surface: r = %+.3f, P = %.2g;  bio05 R2 on lat/lon = %.3f",
      mm, ct$estimate, ct$p.value,
      summary(lm(as.numeric(e$bio05[k]) ~ e$lat[k] + e$lon[k]))$r.squared)
}

# ---- 2. Delta against each predictor, with and without controls ----
TARGETS <- list(
  list(loc = "chr1", pn = "FL1_private",     lab = "CpL1a QUI"),
  list(loc = "chr2_nonFL1_deficit", pn = "non_FL1_deficit", lab = "CpL2a non-QUI"),
  list(loc = "chr3", pn = "CAL1_private",    lab = "CpL3a MOL")
)

for (T in TARGETS) {
  dd <- d[locus == T$loc & panel == T$pn]
  setkey(dd, pool); dd <- dd[e$site]
  y <- dd$delta
  say("")
  say("######## %s ########", T$lab)
  say("%-12s %11s %11s %11s %11s", "predictor", "P plain", "P +metro", "P +lat/lon",
      "adjR2 plain")
  for (v in COVS) {
    x <- as.numeric(e[[v]])
    f1 <- lm(y ~ x)
    f2 <- lm(y ~ x + e$metro)
    f3 <- lm(y ~ x + e$lat + e$lon)
    say("%-12s %11.4g %11.4g %11.4g %11.3f", v,
        summary(f1)$coefficients["x", 4],
        summary(f2)$coefficients["x", 4],
        summary(f3)$coefficients["x", 4],
        summary(f1)$adj.r.squared)
  }
  # bio05 fitted separately within each metro
  say("  --- bio05 within each metro separately ---")
  for (mm in unique(e$metro)) {
    k <- e$metro == mm
    if (sum(k) < 5) { say("    %-11s n=%d, too few", mm, sum(k)); next }
    f <- lm(y[k] ~ as.numeric(e$bio05[k])); s <- summary(f)
    say("    %-11s n=%2d  slope %+.4f  P = %.3g  adjR2 %6.3f  (bio05 span %.2f C)",
        mm, sum(k), coef(f)[2], s$coefficients[2, 4], s$adj.r.squared,
        diff(range(e$bio05[k])))
  }
  kw <- e$metro == "Washington"
  cf <- coef(summary(lm(y[kw] ~ as.numeric(e$bio05[kw]) + as.numeric(e$pis[kw]))))
  say("    Washington, bio05 + impervious surface: P bio05 = %.3g, P pis = %.3g",
      cf[2, 4], cf[3, 4])
}
say("")
say("Reading: 'P +metro' removes the between-metro contrast. If bio05 survives")
say("there, the temperature claim is recoverable. If it collapses while impervious")
say("surface holds, the temperature signal is the Baltimore-vs-Washington contrast.")

# ---- 3. both gradients together, with a metropolitan-area term ----
say("")
say("=== Delta ~ bio05 + pis + metro (both gradients and a metro term) ===")
say("%-34s %10s %10s", "locus / panel", "P bio05", "P pis")
for (lp in unique(paste(d$locus, d$panel))) {
  dd <- d[paste(locus, panel) == lp]
  y  <- dd$delta[match(e$site, dd$pool)]
  cf <- coef(summary(lm(y ~ as.numeric(e$bio05) + as.numeric(e$pis) + e$metro)))
  say("%-34s %10.4f %10.4f", lp, cf[2, 4], cf[3, 4])
}
