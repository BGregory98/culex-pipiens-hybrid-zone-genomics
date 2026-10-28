# =============================================================================
# R/07_multi_locus/25_background_ancestry_vs_winter.R
# =============================================================================
# Whether genome-wide donor ancestry follows winter climate. Form molestus does
# not diapause, and in the Old World it breeds aboveground only where winters
# are mild, so relaxed selection for diapause at warmer sites could raise
# molestus-like ancestry there. The temperature gradient used elsewhere is
# maximum summer temperature (bio05); this script asks whether it also indexes
# winter temperature across the 24 sites, and whether genome-wide alpha for
# each donor follows the minimum temperature of the coldest month (bio06) or the
# mean temperature of the coldest quarter (bio11), alone and with impervious
# surface. Winter layers are CHELSA V2.1 1981-2010 means within the same 500 m
# radius buffer as every other site predictor.
#
# Manuscript: Discussion (molestus-like ancestry, diapause and winter climate)
#
# Inputs:  results/multilocus/genome_wide_panel_vs_env_per_pool.csv  (22_background_ancestry_vs_environment.R)
#          data/site_coords_with_spatial_data_allvars.csv           (R/00_setup/01_site_predictors.R)
# Outputs: results/multilocus/Table_genome_wide_panel_vs_winter.csv
#            (one row per donor x model: focal winter term, slope, P,
#             adjusted R2; plus the correlations of bio05 and impervious
#             surface with each winter variable)
#
# Runtime: seconds
# =============================================================================

source("R/helpers/config.R")
suppressMessages(library(data.table))

PER_POOL <- "results/multilocus/genome_wide_panel_vs_env_per_pool.csv"
ALLVARS  <- "data/site_coords_with_spatial_data_allvars.csv"
need_files(PER_POOL, ALLVARS)
OUT <- "results/multilocus/Table_genome_wide_panel_vs_winter.csv"

pp  <- fread(PER_POOL)
env <- fread(ALLVARS)
d <- merge(pp[, .(pool, donor, alpha_genome, bio05, pis)],
           env[, .(pool = site, bio06, bio11, bio05_allvars = bio05)], by = "pool")
stopifnot(nrow(d) == 72, isTRUE(all.equal(d$bio05, d$bio05_allvars, tolerance = 1e-6)))

# ---- does the summer gradient index winter temperature? ----------------------
sites <- unique(d[, .(pool, bio05, pis, bio06, bio11)])
stopifnot(nrow(sites) == 24)
cors <- rbindlist(lapply(c("bio06", "bio11"), function(w) rbindlist(lapply(c("bio05", "pis"), function(x) {
  ct <- cor.test(sites[[x]], sites[[w]])
  data.table(donor = "sites", model = sprintf("r(%s, %s)", x, w), focal_term = w,
             slope = unname(ct$estimate), p = ct$p.value, adj_r2 = NA_real_, n = 24L)
}))))

# ---- genome-wide alpha against winter temperature ------------------------------
fits <- rbindlist(lapply(c("QUI", "PIP", "MOL"), function(dn) {
  dd <- d[donor == dn]
  rbindlist(lapply(c("bio06", "bio11"), function(w) rbindlist(lapply(c("", " + pis"), function(extra) {
    f <- as.formula(paste("alpha_genome ~", w, extra))
    s <- summary(lm(f, data = dd))
    data.table(donor = dn, model = paste0("alpha ~ ", w, extra), focal_term = w,
               slope = coef(s)[w, 1], p = coef(s)[w, 4], adj_r2 = s$adj.r.squared, n = nrow(dd))
  }))))
}))
out <- rbind(cors, fits)
fwrite(out, OUT)

# ---- console report ----------------------------------------------------
say("Across the 24 sites:")
for (i in seq_len(nrow(cors)))
  say("  %-16s r = %+.2f  (P = %.2g)", cors$model[i], cors$slope[i], cors$p[i])
say("\nGenome-wide alpha against winter temperature (24 field pools):")
say("%-5s %-22s %11s %9s %7s", "donor", "model", "slope", "P", "adjR2")
for (i in seq_len(nrow(fits))) {
  r <- fits[i]
  say("%-5s %-22s %+11.4g %9.2g %7.3f", r$donor, r$model, r$slope, r$p, r$adj_r2)
}
say("\nwrote %s", OUT)
