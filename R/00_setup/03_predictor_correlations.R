# =============================================================================
# R/00_setup/03_predictor_correlations.R
# =============================================================================
# Pairwise correlations among the eight site predictors across the 24 field
# sites (Pearson and Spearman, with 95% CIs and BH-adjusted P over the 28
# pairs), the bio05 x impervious-surface correlation, and the mean absolute
# correlation within the temperature gradient, within the urbanization
# gradient and between them.
#
# Manuscript: Results "Mosquito collection and environmental variation across
#             the corridor" (bio05 vs impervious surface, r = 0.26);
#             Supplementary Methods "Environmental predictors" (bio05-bio07
#             r = 0.84, impervious surface-tree canopy r = -0.90)
#
# Inputs:  data/site_coords_with_spatial_data.csv  (no script: the eight retained
#          predictors of data/site_coords_with_spatial_data_allvars.csv,
#          R/00_setup/01_site_predictors.R, at 10 significant digits)
# Outputs: none (printed)
#
# Runtime: seconds.
# =============================================================================

source("R/helpers/config.R")

PRED_FILE <- 'data/site_coords_with_spatial_data.csv'
need_files(PRED_FILE)

data <- read.csv(PRED_FILE)

## 1. Predictor set ----------------------------------------------------
## Drop identifiers/coords (site, lat, lon). Keep env/landscape covariates.
##   temperature  = bio05 (max temp warmest mo.), bio07 (temp annual range)
##   urbanization = pis (% impervious), tree_canopy (inverse proxy), pd (pop. density)
## elevation, wtd, bio18 retained as context.
env_vars  <- c("elevation","pis","tree_canopy","pd","wtd","bio05","bio07","bio18")
temp_axis <- c("bio05","bio07")
urb_axis  <- c("pis","tree_canopy","pd")

X <- data[, env_vars]
stopifnot(all(vapply(X, is.numeric, logical(1))))
cat("n sites =", nrow(X), "(corridor field sites)\n")

## 2. Correlation matrices --------------------------------------------
## Pearson = linear; Spearman = rank-based, robust to PIS's strong right
## skew / 0-bound. Report both.
R_p <- cor(X, method = "pearson")
R_s <- cor(X, method = "spearman")
cat("\nPearson:\n");  print(round(R_p, 2))
cat("\nSpearman:\n"); print(round(R_s, 2))

## 3. Pairwise r + p + 95% CI (base R; FDR-adjusted over the 28 pairs) --
pair_cor <- function(df, method = "pearson") {
  v <- names(df); rows <- list()
  for (i in 1:(length(v)-1)) for (j in (i+1):length(v)) {
    ct <- suppressWarnings(cor.test(df[[v[i]]], df[[v[j]]], method = method))
    rows[[length(rows)+1]] <- data.frame(
      var1 = v[i], var2 = v[j],
      r  = unname(ct$estimate), p = ct$p.value,
      lo95 = if (!is.null(ct$conf.int)) ct$conf.int[1] else NA_real_,
      hi95 = if (!is.null(ct$conf.int)) ct$conf.int[2] else NA_real_)
  }
  res <- do.call(rbind, rows)
  res$p_BH <- p.adjust(res$p, "BH")
  res[order(abs(res$r), decreasing = TRUE), ]
}
pc <- pair_cor(X, "pearson")
sp <- pair_cor(X, "spearman")            # base cor.test gives no CI for Spearman
cat("\nAll pairs (Pearson, |r|-sorted):\n"); print(pc, row.names = FALSE, digits = 3)

## 4. bio05 x PIS ------------------------------------------------------
gp <- function(tab,a,b) tab[(tab$var1==a&tab$var2==b)|(tab$var1==b&tab$var2==a), ]
cat("\n=== bio05 x pis ===\n")
cat("Pearson :\n");  print(gp(pc,"bio05","pis"), row.names = FALSE, digits = 3)
cat("Spearman:\n");  print(gp(sp,"bio05","pis"), row.names = FALSE, digits = 3)

## 5. Within- and between-gradient summary -----------------------------
## Separable gradients => strong WITHIN-gradient correlation, ~0 BETWEEN.
wt <- R_s[temp_axis, temp_axis]; wt <- mean(abs(wt[upper.tri(wt)]))
wu <- R_s[urb_axis,  urb_axis ]; wu <- mean(abs(wu[upper.tri(wu)]))
bw <- mean(abs(R_s[temp_axis, urb_axis]))
cat(sprintf("\nSpearman mean|r|:  within-temperature %.2f | within-urbanization %.2f | between-gradients %.2f\n",
            wt, wu, bw))
