# =============================================================================
# R/02_qc/01_mapping_rates.R
# =============================================================================
# Tests whether the per-pool mapping rate (fraction of sequenced reads retained
# after MAPQ >= 10, proper-pair and duplicate filtering) covaries with maximum
# summer temperature (bio05) across the 24 field pools: linear regression of
# mapping rate on bio05, plus the range across pools.
#
# Manuscript: Supplementary Methods "Callable-position read-depth analysis"
#             (mapping rates do not covary with bio05, beta = -2.1e-3, P = 0.59)
#
# Inputs:  data/mapping_rates_per_pool.csv          (per-pool read counts from the
#                                                    cluster alignment; no script)
#          data/site_coords_with_spatial_data.csv   (no script: the eight retained
#                                                    predictors of ..._allvars.csv,
#                                                    R/00_setup/01_site_predictors.R)
# Outputs: none (printed)
#
# Runtime: seconds.
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/sample_metadata.R")

library(tidyverse)

MAPPING_FILE <- 'data/mapping_rates_per_pool.csv'
ENV_FILE     <- 'data/site_coords_with_spatial_data.csv'
need_files(MAPPING_FILE, ENV_FILE)

mapping_stats <- read.csv(MAPPING_FILE)

env_data <- read.csv(ENV_FILE) %>%
  dplyr::rename(pool = site)

plot_dat <- mapping_stats %>%
  dplyr::filter(pool %in% field_names) %>%
  dplyr::left_join(env_data %>% dplyr::select(pool, bio05), by = 'pool')

# Regression of mapping rate on bio05
fit_mapping <- lm(mapping_rate ~ bio05, data = plot_dat)
print(summary(fit_mapping))

# Across-pool variation
cat(sprintf("Mapping rate range: %.3f - %.3f (max %.1f%% across-pool variation)\n",
            min(plot_dat$mapping_rate),
            max(plot_dat$mapping_rate),
            100 * (max(plot_dat$mapping_rate) - min(plot_dat$mapping_rate))))
