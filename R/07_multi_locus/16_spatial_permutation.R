# =============================================================================
# R/07_multi_locus/16_spatial_permutation.R
# =============================================================================
# Spatial autocorrelation of the two predictors, and a permutation test of the
# private-panel ancestry-environment correlations that tries to preserve it.
# Moran's I is computed on row-standardized inverse great-circle-distance
# weights. bio05 is largely predicted by latitude and longitude, so a plain
# shuffle across the 24 sites destroys its spatial structure and makes the null
# too easy (anti-conservative). Shuffles are therefore retained only if their
# Moran's I lies within 0.15 of the observed predictor's; the panel alpha
# correlations are tested against the retained shuffles and, for comparison,
# against unconstrained ones.
#
# Manuscript: Supplementary Methods "Geographic and metropolitan-area controls"
#             (observed Moran's I of bio05, and its range across 2,000
#             unconstrained permutations)
#
# Inputs:  results/multilocus/Table_private_panel_alpha_per_pool.csv  (R/07_multi_locus/12_private_panel_clines.R)
#          data/site_coords_with_spatial_data_allvars.csv             (R/00_setup/01_site_predictors.R)
# Outputs: none (console only)
#
# Runtime: a few minutes (up to 1,000,000 shuffles per predictor)
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/geo.R")
suppressMessages(library(data.table))

need_files("results/multilocus/Table_private_panel_alpha_per_pool.csv",
           "data/site_coords_with_spatial_data_allvars.csv")

set.seed(42)

N_KEEP  <- 2000    # retained permutations per predictor
MORAN_TOL <- 0.15  # |I_perm - I_obs| tolerance, on Moran's I

alpha <- fread("results/multilocus/Table_private_panel_alpha_per_pool.csv", showProgress = FALSE)
env   <- fread("data/site_coords_with_spatial_data_allvars.csv", showProgress = FALSE)

sites <- sort(unique(alpha$pool))
setkey(env, site); e <- env[sites]
stopifnot(nrow(e) == length(sites))

# Inverse-distance spatial weights, row-standardised
n <- nrow(e)
D <- outer(seq_len(n), seq_len(n), Vectorize(function(i, j)
  if (i == j) 0 else hav(e$lat[i], e$lon[i], e$lat[j], e$lon[j])))
W <- ifelse(D > 0, 1 / D, 0)
W <- W / rowSums(W)

moran <- function(x) {
  z <- x - mean(x)
  (n / sum(W)) * sum(W * outer(z, z)) / sum(z^2)
}

say("sites: %d   mean pairwise distance: %.1f km", n, mean(D[upper.tri(D)]))

for (predictor in c("bio05", "pis")) {
  x_obs <- as.numeric(e[[predictor]])
  I_obs <- moran(x_obs)
  say("")
  say("=== %s: observed Moran's I = %+.3f ===", predictor, I_obs)

  # rejection-sample permutations that retain comparable spatial structure
  keep <- vector("list", N_KEEP); k <- 0; tries <- 0
  naive_I <- numeric(0)
  while (k < N_KEEP && tries < N_KEEP * 500) {
    tries <- tries + 1
    xp <- sample(x_obs)
    Ip <- moran(xp)
    if (tries <= 2000) naive_I <- c(naive_I, Ip)
    if (abs(Ip - I_obs) <= MORAN_TOL) { k <- k + 1; keep[[k]] <- xp }
  }
  keep <- keep[seq_len(k)]
  say("  retained %d of %d shuffles (%.2f%%) with |I - I_obs| <= %.2f",
      k, tries, 100 * k / tries, MORAN_TOL)
  say("  unconstrained shuffles give Moran's I = %+.3f (mean), %+.3f to %+.3f",
      mean(naive_I), min(naive_I), max(naive_I))
  say("  -> a naive shuffle would test against near-zero spatial structure,")
  say("     which is why it is anti-conservative here.")
  say("  retained permutations: Moran's I mean %+.3f",
      mean(sapply(keep, moran)))

  for (loc in unique(alpha$locus)) {
    for (pn in unique(alpha$panel[alpha$locus == loc])) {
      d <- alpha[locus == loc & panel == pn]
      setkey(d, pool); d <- d[sites]
      if (sd(d$alpha) == 0) next
      # only the predictor tested at each locus
      wanted <- (loc %in% c("CpL1a", "CpL2a") && predictor == "bio05") ||
                (loc == "CpL3a" && predictor == "pis")
      if (!wanted) next
      r_obs <- cor(d$alpha, x_obs)
      r_null <- sapply(keep, function(xp) cor(d$alpha, xp))
      p_abs <- mean(abs(r_null) >= abs(r_obs))
      p_dir <- if (r_obs >= 0) mean(r_null >= r_obs) else mean(r_null <= r_obs)
      r_naive <- replicate(2000, cor(d$alpha, sample(x_obs)))
      say("  %-6s %-12s r = %+.3f | spatial-constrained P(|r|) = %.4f  dir = %.4f | naive-shuffle P(|r|) = %.4f",
          loc, pn, r_obs, p_abs, p_dir, mean(abs(r_naive) >= abs(r_obs)))
    }
  }
}
say("")
say("done.")
