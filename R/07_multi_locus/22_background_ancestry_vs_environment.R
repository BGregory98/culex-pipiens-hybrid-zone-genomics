# =============================================================================
# R/07_multi_locus/22_background_ancestry_vs_environment.R
# =============================================================================
# Whether the GENOME-WIDE frequency of each donor-private panel (alpha_genome,
# the background against which every locus-level Delta is measured) tracks the
# environment across the 24 field pools. alpha_genome is regressed on bio05,
# impervious surface, both, bio05 with a metropolitan-area term, and bio05 with
# pool heterozygosity, so the background and the loci are on one scale. The
# impervious-surface association is then refitted against position (coordinates,
# distance to central Washington), within each metro area, and without the four
# central-Washington pools.
#
# Manuscript: Results "Genome-wide admixture is structured by both temperature
#             and urbanization" (molestus-like ancestry against impervious
#             surface, with its position and metro-area controls); Supplementary
#             Results "Temperature is collinear with position along the corridor"
#             (genome-wide QUI-private ancestry against bio05 and MOL-private
#             ancestry against impervious surface) and "Genome-wide molestus-like
#             ancestry follows impervious surface independently of position";
#             Discussion (urbanization)
#
# Inputs:  results/multilocus/per_pool_ancestry_per_pool.csv  (R/07_multi_locus/01_per_pool_ancestry.R)
#          data/site_coords_with_spatial_data.csv             (R/00_setup/01_site_predictors.R)
#          results/multilocus/Table_per_pool_diversity.csv    (R/04_pop_structure/05_pool_diversity.R)
# Outputs: results/multilocus/Table_genome_wide_panel_vs_env.csv
#          results/multilocus/genome_wide_panel_vs_env_per_pool.csv
#          results/multilocus/Table_genome_wide_panel_vs_pis_robustness.csv
#
# Runtime: seconds
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/geo.R")
suppressMessages(library(data.table))

need_files("results/multilocus/per_pool_ancestry_per_pool.csv",
           "data/site_coords_with_spatial_data.csv",
           "results/multilocus/Table_per_pool_diversity.csv")

# ---- load ---------------------------------------------------------------
pp  <- fread("results/multilocus/per_pool_ancestry_per_pool.csv")
env <- fread("data/site_coords_with_spatial_data.csv")
div <- fread("results/multilocus/Table_per_pool_diversity.csv")

# genome-wide alpha is identical across loci for a given panel; take one
gw <- pp[panel %in% c("FL1_private", "ME1_private", "CAL1_private"),
         .(pool, panel, alpha_genome)]
gw <- unique(gw)
stopifnot(nrow(gw) == 72)
gw[, donor := fcase(panel == "FL1_private", "QUI",
                    panel == "ME1_private", "PIP",
                    panel == "CAL1_private", "MOL")]

# ---- metro assignment (nearer city centre, as in 18_within_metro_models.R) ----
env[, d_bal := hav(lat, lon, BAL[1], BAL[2])]
env[, d_dc  := hav(lat, lon, DC[1],  DC[2])]
env[, metro := ifelse(d_bal < d_dc, "Baltimore", "Washington")]

d <- merge(gw, env[, .(pool = site, bio05, pis, tree_canopy, lat, lon, metro)],
           by = "pool")
d <- merge(d, div[, .(pool, exp_het)], by = "pool")
stopifnot(nrow(d) == 72)
fwrite(d[order(donor, bio05)], "results/multilocus/genome_wide_panel_vs_env_per_pool.csv")

# ---- models ---------------------------------------------------------------
fit_row <- function(dd, donor, label, formula) {
  m  <- lm(formula, data = dd)
  s  <- summary(m)
  co <- coef(s)
  # report the first non-intercept term as the focal predictor
  focal <- rownames(co)[2]
  data.table(donor = donor, model = label, focal_term = focal,
             slope = co[2, 1], se = co[2, 2], t = co[2, 3], p = co[2, 4],
             adj_r2 = s$adj.r.squared, n = nrow(dd),
             other_terms = paste(rownames(co)[-(1:2)], collapse = ";"),
             other_p = paste(signif(co[-(1:2), 4], 3), collapse = ";"))
}

out <- rbindlist(lapply(c("QUI", "PIP", "MOL"), function(dn) {
  dd <- d[donor == dn]
  rbindlist(list(
    fit_row(dd, dn, "alpha ~ bio05",                 alpha_genome ~ bio05),
    fit_row(dd, dn, "alpha ~ pis",                   alpha_genome ~ pis),
    fit_row(dd, dn, "alpha ~ bio05 + pis",           alpha_genome ~ bio05 + pis),
    fit_row(dd, dn, "alpha ~ pis + bio05",           alpha_genome ~ pis + bio05),
    fit_row(dd, dn, "alpha ~ bio05 + metro",         alpha_genome ~ bio05 + metro),
    fit_row(dd, dn, "alpha ~ bio05 + exp_het",       alpha_genome ~ bio05 + exp_het),
    fit_row(dd, dn, "alpha ~ bio05 + pis + exp_het", alpha_genome ~ bio05 + pis + exp_het)
  ))
}))
fwrite(out, "results/multilocus/Table_genome_wide_panel_vs_env.csv")

# ---- console report ----------------------------------------------------
say("\nGenome-wide donor-private panel alpha vs environment (24 field pools)")
say("%-5s %-28s %-8s %10s %9s %8s %7s  %s", "donor", "model", "focal",
    "slope", "se", "p", "adjR2", "other terms (p)")
for (i in seq_len(nrow(out))) {
  r <- out[i]
  say("%-5s %-28s %-8s %10.4f %9.4f %8.4f %7.3f  %s (%s)", r$donor, r$model,
      r$focal_term, r$slope, r$se, r$p, r$adj_r2, r$other_terms, r$other_p)
}

# range of alpha across the corridor per donor, to size the effect
say("\nAlpha range and predicted change over observed bio05 span:")
bspan <- diff(range(d$bio05))
for (dn in c("QUI", "PIP", "MOL")) {
  dd <- d[donor == dn]
  m  <- lm(alpha_genome ~ bio05, dd)
  say("  %s: observed alpha %.3f-%.3f (mean %.3f); slope %+.4f per degC x %.2f degC span = %+.4f",
      dn, min(dd$alpha_genome), max(dd$alpha_genome), mean(dd$alpha_genome),
      coef(m)[2], bspan, coef(m)[2] * bspan)
}

# cool-six vs warm-six means, matching how the locus results are reported
say("\nSix coolest vs six warmest sites (mean alpha):")
setorder(d, bio05)
for (dn in c("QUI", "PIP", "MOL")) {
  dd <- d[donor == dn]
  say("  %s: coolest six %.3f, warmest six %.3f",
      dn, mean(head(dd$alpha_genome, 6)), mean(tail(dd$alpha_genome, 6)))
}

# ---- urbanization: is the impervious-surface association a position effect? ----
# Impervious surface has no spatial trend across the corridor, unlike bio05, so
# a genome-wide association with it cannot come from proximity to a donor's
# source region. Each donor's alpha ~ pis is refitted with the coordinates, with
# distance to central Washington, with bio05 and metro area, within each metro
# area, and without the four pools within 3 km of central Washington (the next
# nearest pool is 11 km away), which sit at the urban end of the gradient.
d <- merge(d, env[, .(pool = site, d_dc)], by = "pool")
CORE_KM <- 3
rob <- rbindlist(lapply(c("QUI", "PIP", "MOL"), function(dn) {
  dd <- d[donor == dn]
  fits <- list(
    "alpha ~ pis"                    = alpha_genome ~ pis,
    "alpha ~ pis + lat + lon"        = alpha_genome ~ pis + lat + lon,
    "alpha ~ pis + dist_dc"          = alpha_genome ~ pis + d_dc,
    "alpha ~ pis + bio05 + metro"    = alpha_genome ~ pis + bio05 + metro)
  r1 <- rbindlist(lapply(names(fits), function(k) {
    s <- summary(lm(fits[[k]], data = dd))
    data.table(donor = dn, test = k, n = nrow(dd), estimate = coef(s)["pis", 1],
               p = coef(s)["pis", 4], adj_r2 = s$adj.r.squared)
  }))
  sub <- list(Baltimore = dd[metro == "Baltimore"], Washington = dd[metro == "Washington"],
              without_dc_core = dd[d_dc > CORE_KM])
  r2 <- rbindlist(lapply(names(sub), function(k) {
    ct <- cor.test(sub[[k]]$alpha_genome, sub[[k]]$pis)
    data.table(donor = dn, test = paste0("r(alpha, pis), ", k), n = nrow(sub[[k]]),
               estimate = unname(ct$estimate), p = ct$p.value, adj_r2 = NA_real_)
  }))
  rbind(r1, r2)
}))
fwrite(rob, "results/multilocus/Table_genome_wide_panel_vs_pis_robustness.csv")
s_pis <- summary(lm(pis ~ lat + lon, data = d[donor == "MOL"]))
say("\nImpervious surface on the coordinates: R2 = %.3f; r with distance to central Washington = %+.2f",
    s_pis$r.squared, cor(d[donor == "MOL"]$pis, d[donor == "MOL"]$d_dc))
say("Genome-wide alpha ~ impervious surface, robustness (the pis term, or r within a subset):")
for (i in seq_len(nrow(rob)))
  say("  %-4s %-34s n=%2d  estimate %+.4g  P = %.3g", rob$donor[i], rob$test[i],
      rob$n[i], rob$estimate[i], rob$p[i])

say("\nwrote results/multilocus/Table_genome_wide_panel_vs_env.csv")
say("wrote results/multilocus/genome_wide_panel_vs_env_per_pool.csv")
say("wrote results/multilocus/Table_genome_wide_panel_vs_pis_robustness.csv")
