# =============================================================================
# R/00_setup/05_table_s1_predictors.R
# =============================================================================
# Table S1 values at full precision. For each of the 24 field sites: its
# metropolitan area (the nearer of the Baltimore and Washington city centres by
# great-circle distance), urban class, coordinates, bio05 and percent impervious
# surface exactly as they enter the models, and the distances to both centres.
# This is the unrounded table that the Table S1 legend says is deposited with
# the code.
#
# Manuscript: Table S1; Supplementary Methods "Geographic and metropolitan-area
#             controls" (assignment rule; 11 Baltimore and 13 Washington sites);
#             read by 06_temperature_vs_impervious.R for Figure S1
#
# Inputs:  data/site_coords_with_spatial_data_allvars.csv   (R/00_setup/01_site_predictors.R)
#          data/sample_categories.csv
# Outputs: results/multilocus/Table_S1_predictors_full_precision.csv
#
# Runtime: seconds
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/sample_metadata.R")
source("R/helpers/geo.R")
suppressMessages(library(data.table))

need_files("data/site_coords_with_spatial_data_allvars.csv")

env <- fread("data/site_coords_with_spatial_data_allvars.csv", showProgress = FALSE)
setkey(env, site); e <- env[field_names]
e$class <- sample_info_all$sample_type_2[match(e$site, sample_info_all$sample)]

e$d_bal <- hav(e$lat, e$lon, BAL[1], BAL[2])
e$d_dc  <- hav(e$lat, e$lon, DC[1],  DC[2])
e$metro <- ifelse(e$d_bal < e$d_dc, "Baltimore", "Washington")
say("metropolitan split: Baltimore %d, Washington %d",
    sum(e$metro == "Baltimore"), sum(e$metro == "Washington"))

# Values are written as stored; nothing is rounded.
dir.create("results/multilocus", showWarnings = FALSE, recursive = TRUE)
fwrite(e[, .(site, metro, class, lat, lon, bio05, pis,
             dist_baltimore_km = d_bal, dist_washington_km = d_dc)],
       "results/multilocus/Table_S1_predictors_full_precision.csv")
say("wrote results/multilocus/Table_S1_predictors_full_precision.csv (%d sites)", nrow(e))
