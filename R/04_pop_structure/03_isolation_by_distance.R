# =============================================================================
# R/04_pop_structure/03_isolation_by_distance.R
# =============================================================================
# Pairwise FST heatmap among the 24 field pools, ordered by maximum summer
# temperature (bio05), and an isolation-by-distance test: one-sided Pearson
# correlation and OLS regression of linearized FST/(1 - FST) (Rousset 1997) on
# log10 geodesic distance over all 276 pairs, with a Mantel permutation test
# (9,999 permutations) as the significance test that allows for the
# non-independence of pairs.
#
# Manuscript: Figure S3; Figure S6; Results "mean pairwise FST 0.004, range
#             0.0003-0.0110" and the isolation-by-distance paragraph
#             (Methods "Population structure, differentiation and isolation by
#             distance")
#
# Inputs:  data/pairwise_fst_all.csv                  (R/04_pop_structure/01_pairwise_fst.R)
#          data/site_coords_with_spatial_data.csv     (R/00_setup/01_site_predictors.R)
#          data/sample_categories.csv                 (sample metadata)
# Outputs: figures/final/supplement/Figure_S3.png     (pairwise FST heatmap)
#          figures/final/supplement/Figure_S6.png     (isolation by distance)
#
# Runtime: < 1 min
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/sample_metadata.R")   # field_names

suppressPackageStartupMessages({
  library(tidyverse)   # provides dplyr, ggplot2, tibble, ...
  library(sf)
  library(viridis)
  library(vegan)
})

set.seed(42)   # Mantel permutations

SITE_COORDS_FILE  <- "data/site_coords_with_spatial_data.csv"
PAIRWISE_FST_FILE <- "data/pairwise_fst_all.csv"
need_files(SITE_COORDS_FILE, PAIRWISE_FST_FILE)

FIG_HEATMAP <- "figures/final/supplement/Figure_S3.png"
FIG_IBD     <- "figures/final/supplement/Figure_S6.png"
dir.create(dirname(FIG_HEATMAP), recursive = TRUE, showWarnings = FALSE)

FST_HEATMAP_LIMIT <- 0.015       # upper limit of FST colour scale

# ---- 1. Field-site coordinates ---------------------------------------------
site_coords <- read.csv(SITE_COORDS_FILE)[, 1:3] %>%
  filter(site %in% field_names) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Sites ordered by bio05 (500 m buffer means, Table S1) for the heatmap axis
site_env    <- read.csv(SITE_COORDS_FILE)
bio05_order <- site_coords$site[order(site_env$bio05[match(site_coords$site, site_env$site)])]

# ---- 2. Load pairwise FST and build the 24x24 field matrix ----------------
pairwise_long <- read.csv(PAIRWISE_FST_FILE)

# Sanity: file should be in {pop1, pop2, fst} long format.
stopifnot(all(c("pop1", "pop2", "fst") %in% colnames(pairwise_long)))

# Keep only field-vs-field pairs.
pairwise_field <- pairwise_long %>%
  filter(pop1 %in% field_names, pop2 %in% field_names)

# Symmetrize: 01_pairwise_fst.R wrote the lower triangle only.
# Build the 24x24 symmetric matrix by populating both (i,j) and (j,i).
pairwise_fst_matrix <- matrix(NA_real_,
                              nrow = length(field_names),
                              ncol = length(field_names),
                              dimnames = list(field_names, field_names))
for (k in seq_len(nrow(pairwise_field))) {
  i <- pairwise_field$pop1[k]; j <- pairwise_field$pop2[k]; v <- pairwise_field$fst[k]
  pairwise_fst_matrix[i, j] <- v
  pairwise_fst_matrix[j, i] <- v
}
diag(pairwise_fst_matrix) <- 0

n_field_pairs <- choose(length(field_names), 2)
if (sum(!is.na(pairwise_fst_matrix)) - length(field_names) != 2 * n_field_pairs) {
  stop("Pairwise FST matrix incomplete after symmetrization — check pairwise_fst_all.csv.")
}

# ---- 3. Figure S3: heatmap of pairwise FST (ordered by bio05) --------------
fst_long <- pairwise_fst_matrix[bio05_order, bio05_order] %>%
  as.table() %>%
  as.data.frame() %>%
  dplyr::rename(pop1 = Var1, pop2 = Var2, fst = Freq) %>%
  filter(match(pop1, bio05_order) > match(pop2, bio05_order)) %>%
  mutate(
    pop1 = factor(pop1, levels = bio05_order),
    pop2 = factor(pop2, levels = bio05_order)
  )

p_heatmap <- ggplot(fst_long, aes(x = pop2, y = pop1, fill = fst)) +
  geom_tile() +
  scale_fill_viridis(option = "B", name = expression(F[ST]),
                     limits = c(0, FST_HEATMAP_LIMIT)) +
  coord_equal() +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title  = element_blank())

ggsave(FIG_HEATMAP, p_heatmap, width = 8, height = 8, units = "in")

# ---- 4. Isolation-by-distance ----------------------------------------------
# Reorder coordinates to match the FST matrix
site_coords_field <- site_coords[match(rownames(pairwise_fst_matrix),
                                       site_coords$site), ]
stopifnot(all(site_coords_field$site == rownames(pairwise_fst_matrix)))

# Pairwise geodesic distances in km (st_distance returns metres on geographic CRS)
dist_km <- units::drop_units(sf::st_distance(site_coords_field)) / 1000
rownames(dist_km) <- colnames(dist_km) <- rownames(pairwise_fst_matrix)

# Lower-triangle long-format table
idx <- which(lower.tri(pairwise_fst_matrix), arr.ind = TRUE)
ibd_df <- tibble(
  pop1     = rownames(pairwise_fst_matrix)[idx[, 1]],
  pop2     = colnames(pairwise_fst_matrix)[idx[, 2]],
  fst      = pairwise_fst_matrix[idx],
  dist_km  = dist_km[idx],
  fst_lin  = fst / (1 - fst),
  log_dist = log10(dist_km)
)
stopifnot(nrow(ibd_df) == choose(nrow(pairwise_fst_matrix), 2))

# One-sided Pearson correlation (IBD predicts r > 0)
ibd_cor <- cor.test(ibd_df$fst_lin, ibd_df$log_dist,
                    method = "pearson", alternative = "greater")

# Linear regression for slope
ibd_lm <- lm(fst_lin ~ log_dist, data = ibd_df)

# Mantel permutation test (accounts for non-independence of pairs)
fst_lin_mat <- pairwise_fst_matrix / (1 - pairwise_fst_matrix)
diag(fst_lin_mat) <- 0
mantel_res <- vegan::mantel(as.dist(fst_lin_mat), as.dist(log10(dist_km)),
                            method = "pearson", permutations = 9999,
                            na.rm  = TRUE)
mantel_r <- mantel_res$statistic
mantel_p <- mantel_res$signif

# Print summary
cat("\n=== Isolation-by-distance: 24 field pools (n =", nrow(ibd_df),
    "pairs) ===\n",
    sprintf("Distance range:           %.1f - %.1f km\n",
            min(ibd_df$dist_km), max(ibd_df$dist_km)),
    sprintf("FST range:                %.4f - %.4f\n",
            min(ibd_df$fst), max(ibd_df$fst)),
    sprintf("Pearson r:                %.4f\n", ibd_cor$estimate),
    sprintf("  one-sided P (>0):       %.4g\n", ibd_cor$p.value),
    sprintf("LM slope:                 %.4g (SE %.4g)\n",
            coef(ibd_lm)[2], summary(ibd_lm)$coefficients[2, 2]),
    sprintf("LM adj. R-squared:        %.4f\n", summary(ibd_lm)$adj.r.squared),
    sprintf("Mantel r (9999 perms):    %.4f (P = %.4g)\n",
            mantel_r, mantel_p),
    "\n",
    sprintf("Manuscript: r = %.3f; one-sided P = %.3g; slope = %.3g\n",
            ibd_cor$estimate, ibd_cor$p.value, coef(ibd_lm)[2]),
    "===================================================================\n\n",
    sep = "")

# ---- 5. Figure S6: IBD plot ------------------------------------------------
p_ibd <- ggplot(ibd_df, aes(x = log_dist, y = fst_lin)) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "grey30", fill = "grey80", linewidth = 0.6) +
  geom_point(alpha = 0.6, size = 1.8) +
  labs(x = expression(log[10] ~ "pairwise geographic distance (km)"),
       y = expression(F[ST] / (1 - F[ST]))) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(FIG_IBD, p_ibd, width = 5, height = 4, units = "in", dpi = 300)
