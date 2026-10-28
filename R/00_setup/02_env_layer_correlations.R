# =============================================================================
# R/00_setup/02_env_layer_correlations.R
# =============================================================================
# Pearson correlations among the 15 candidate environmental layers across all
# raster cells of the study area, used to prune the exploratory set to the
# eight predictors of the genome scan. Lists every pair with |r| > 0.8 and the
# retained predictors.
#
# Manuscript: Supplementary Methods "Environmental predictors" (pruning for
#             redundancy; bio06 and bio11 dropped for > 0.9 correlation with
#             elevation)
#
# Inputs:  data/local_env_layers_2.tif   (R/00_setup/01_site_predictors.R)
# Outputs: none (printed)
#
# Runtime: seconds.
# =============================================================================

source("R/helpers/config.R")

library(terra)
library(tibble)
library(dplyr)

RASTER_FILE <- 'data/local_env_layers_2.tif'
need_files(RASTER_FILE)

spatial_data <- rast(RASTER_FILE)

df <- as.data.frame(spatial_data, na.rm = TRUE)

cor_mat <- cor(df, use = "pairwise.complete.obs", method = "pearson")

# Pairs with |r| > 0.8
high_corr <- which(abs(cor_mat) > 0.8 & abs(cor_mat) < 1, arr.ind = TRUE)

high_corr_tbl <- as_tibble(high_corr) %>%
  mutate(
    var1 = rownames(cor_mat)[row],
    var2 = colnames(cor_mat)[col],
    corr = cor_mat[cbind(row, col)]
  ) %>%
  filter(var1 < var2)  # remove duplicates

cat('\nLayer pairs with |r| > 0.8 across raster cells:\n')
print(as.data.frame(high_corr_tbl %>% dplyr::select(var1, var2, corr) %>%
                      arrange(desc(abs(corr)))), row.names = FALSE, digits = 3)

variables_to_remove = c('bio06', 'bio11', 'fcf', 'gdd10',  # all 90%+ correlated with elevation
                        'ndvi',            # 85-90% correlated with tree_canopy and pis
                        'hd',              # 90% correlated with pd
                        'bio04'            # 74% correlated with bio07 and bio18
                      )

# Bio05 = highest monthly mean of daily maximum temperatures; peak thermal conditions
# Bio07 = annual temperature range (warmest minus coldest month)
# Bio18 = mean monthly precipitation of the warmest quarter
kept_variables = names(spatial_data)[!names(spatial_data) %in% variables_to_remove]
say("\nRetained predictors (%d): %s", length(kept_variables), paste(kept_variables, collapse = ", "))
