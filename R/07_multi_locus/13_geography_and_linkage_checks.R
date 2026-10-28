# =============================================================================
# R/07_multi_locus/13_geography_and_linkage_checks.R
# =============================================================================
# Two checks on tables already in the repository.
#   1. How far maximum summer temperature (bio05) is explained by position along
#      the corridor: bio05 ~ latitude + longitude, its residual SD, and bio05's
#      correlation with longitude, latitude and great-circle distance to the QUI
#      reference source (Vero Beach, Table S9); the warmest sites by urban class.
#   2. Linkage amplitude against distance: for each locus and pair set, median
#      r2 at the shortest and longest distance bins, and the flank-pair curve
#      against the same-chromosome baseline.
#
# Manuscript: Supplementary Results "Temperature is collinear with position along
#             the corridor" (R2 = 0.858, residual and raw SD, r with longitude
#             and with distance to the QUI source); Results "Mosquito collection
#             and environmental variation across the corridor" (r = 0.26);
#             Results CpL1a section (the four warmest urban Washington sites);
#             Results CpL2a section and Supplementary Results "Linkage
#             architecture at the three candidate loci" (flank r2 by distance)
#
# Inputs:  data/site_coords_with_spatial_data_allvars.csv   (R/00_setup/01_site_predictors.R)
#          data/sample_categories.csv
#          results/LD/Table_multi_locus_LD_decay.csv        (R/07_multi_locus/04_linkage_decay.R)
# Outputs: none (console only)
#
# Runtime: seconds
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/sample_metadata.R")
source("R/helpers/geo.R")

need_files("data/site_coords_with_spatial_data_allvars.csv",
           "results/LD/Table_multi_locus_LD_decay.csv")

# -----------------------------------------------------------------------------
# 1. bio05 against geography
# -----------------------------------------------------------------------------
env <- read.csv("data/site_coords_with_spatial_data_allvars.csv", stringsAsFactors = FALSE)
e <- env[env$site %in% field_names, ]
say("field sites: %d", nrow(e))

say("")
say("=== bio05 against geography ===")
fit <- lm(bio05 ~ lat + lon, data = e); s <- summary(fit)
say("bio05 ~ lat + lon      : R2 = %.3f (adj %.3f), P = %.3g",
    s$r.squared, s$adj.r.squared, pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE))
say("raw SD(bio05)          : %.3f C", sd(e$bio05))
say("residual SD after geog : %.3f C", sd(resid(fit)))
say("r(bio05, longitude)    : %+.3f", cor(e$bio05, e$lon))
say("r(bio05, latitude)     : %+.3f", cor(e$bio05, e$lat))

# Great-circle distance to the Cx. quinquefasciatus source (FL1):
# Vero Beach, FL, the QUI collection site in Table S9.
VB <- c(lat = 27.568669, lon = -80.417803)
e$dist_qui <- hav(e$lat, e$lon, VB["lat"], VB["lon"])
say("r(bio05, distance to QUI source): %+.3f   [Vero Beach, FL]", cor(e$bio05, e$dist_qui))
say("r(impervious surface, bio05)    : %+.3f", cor(e$pis, e$bio05))
say("pairwise site distance range    : %.1f to %.1f km",
    min(dist(cbind(e$lat*111, e$lon*111*cos(mean(e$lat)*pi/180)))),
    max(dist(cbind(e$lat*111, e$lon*111*cos(mean(e$lat)*pi/180)))))

say("")
say("=== the warmest sites and their urban class ===")
o <- e[order(-e$bio05), c("site", "bio05")]
o$class <- sample_info_all$sample_type_2[match(o$site, sample_info_all$sample)]
print(head(o, 7), row.names = FALSE)
say("field_u = urban, field_s = suburban, field_r = rural (data/sample_categories.csv)")

# -----------------------------------------------------------------------------
# 2. Linkage amplitude against distance
# -----------------------------------------------------------------------------
ld <- read.csv("results/LD/Table_multi_locus_LD_decay.csv", stringsAsFactors = FALSE)
say("")
say("=== linkage amplitude against distance: LD sets available ===")
print(table(ld$locus, ld$set))

# order distance bins by their geometric midpoint
bin_mid <- function(b) {
  b <- gsub("kb", "e3", b); b <- gsub("Mb", "e6", b)
  p <- strsplit(b, "-")[[1]]
  lo <- as.numeric(sub("e3$", "", sub("e6$", "", p[1])))
  # scale of the lower bound follows the upper bound's unit
  unit <- if (grepl("e6", p[2])) 1e6 else 1e3
  hi <- as.numeric(sub("e3$", "", sub("e6$", "", p[2]))) * unit
  lo <- lo * unit
  sqrt(max(lo, 1) * hi)
}
ld$mid <- sapply(ld$distance_bin, bin_mid)

for (L in c("chr1", "chr2", "chr3")) {
  d <- ld[ld$locus == L, ]
  say("")
  say("---- %s ----", L)
  for (st in sort(unique(d$set))) {
    dd <- d[d$set == st, ]
    dd <- dd[order(dd$mid), ]
    say("  %-34s median r2 %.3f -> %.3f across %s .. %s",
        st, dd$median_r2[1], dd$median_r2[nrow(dd)],
        dd$distance_bin[1], dd$distance_bin[nrow(dd)])
  }
  # flank-pair (set C) curve in full, plus the baseline it is compared against
  fl <- d[grepl("^C", d$set), ]; fl <- fl[order(fl$mid), ]
  bl <- d[grepl("^E", d$set), ]
  if (nrow(fl)) {
    say("  flank pairs by distance:")
    for (i in seq_len(nrow(fl))) say("     %-12s %.3f", fl$distance_bin[i], fl$median_r2[i])
    if (nrow(bl)) say("  same-chromosome baseline: %.3f", median(bl$median_r2))
  }
}
say("")
say("done.")
