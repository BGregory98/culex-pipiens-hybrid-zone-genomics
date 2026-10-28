# =============================================================================
# R/00_setup/01_site_predictors.R
# =============================================================================
# Builds the environmental raster stack for the study area and extracts the
# site predictors. Seven landscape layers (elevation, impervious surface, tree
# canopy, NDVI, population density, water-table depth, housing density) are
# joined by eight CHELSA V2.1 climate layers cropped to the site bounding box
# (+/- 0.2 degrees) and resampled to the landscape grid. Each layer is then
# projected to UTM 18N and averaged within a 500 m radius (1 km diameter)
# buffer around each of the 24 field sites.
#
# Manuscript: Methods "Population structure, differentiation and isolation by
#             distance" and Supplementary Methods "Environmental predictors"
#             (predictors extracted within a 1 km diameter buffer); Table S1
#
# Inputs:  data/dc_env_layers.tif                          (landscape stack; no script)
#          data/site_coords_with_buffered_spatial_data.csv (site coordinates; no script)
#          deposit: data/chelsa_layers/CHELSA_<var>_1981-2010_V.2.1.tif, <var> =
#                   bio04 bio05 bio06 bio07 bio11 bio18 fcf gdd10
#                   (public CHELSA V2.1 climatologies 1981-2010; https://chelsa-climate.org)
# Outputs: data/local_env_layers_2.tif                     (15-layer stack)
#          data/site_coords_with_spatial_data_allvars.csv  (24 sites x 15 predictors)
#          data/site_coords_with_spatial_data.csv          (the 8 retained predictors,
#                                                           at most 11 characters per value;
#                                                           also the BayPass covariate table)
#
# Runtime: a few minutes (reads the global CHELSA rasters).
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/sample_metadata.R")

library(sf)
library(terra)
library(raster)
library(dplyr)

LANDSCAPE_FILE <- 'data/dc_env_layers.tif'
SITES_FILE     <- 'data/site_coords_with_buffered_spatial_data.csv'
need_files(LANDSCAPE_FILE, SITES_FILE)

CHELSA_VARS  <- c('bio04', 'bio05', 'bio06', 'bio07', 'bio11', 'bio18', 'fcf', 'gdd10')
chelsa_files <- vapply(sprintf('data/chelsa_layers/CHELSA_%s_1981-2010_V.2.1.tif', CHELSA_VARS),
                       deposit_file, character(1), USE.NAMES = FALSE)

OUT_STACK <- 'data/local_env_layers_2.tif'
OUT_TABLE <- 'data/site_coords_with_spatial_data_allvars.csv'
OUT_KEPT  <- 'data/site_coords_with_spatial_data.csv'

# Predictors kept after the redundancy pruning (02_env_layer_correlations.R).
KEPT <- c('elevation', 'pis', 'tree_canopy', 'pd', 'wtd', 'bio05', 'bio07', 'bio18')

### MAP ###
# Landscape stack without layers 3 (land_cover) and 9 (dcc).
spatial_data <- rast(LANDSCAPE_FILE)[[-c(3,9)]]

# Site coordinates
site_coords = read.csv(SITES_FILE)[,1:3]
site_coords = site_coords %>% arrange(site)
site_coords = site_coords %>% filter(site %in% field_names)

# sf object from lat/lon (coords = c(longitude, latitude))
sites_sf <- st_as_sf(site_coords, coords = c("lon", "lat"), crs = 'EPSG:4326')

# Bounding box of the sites +/- 0.2 degrees (WGS84)
borders <- c(xmin = min(site_coords$lon) - 0.2, xmax = max(site_coords$lon) + 0.2,
             ymin = min(site_coords$lat) - 0.2, ymax = max(site_coords$lat) + 0.2)
bbox_wgs84 <- st_as_sfc(st_bbox(borders, crs = crs(sites_sf)))
extent_smaller <- st_sf(geometry = bbox_wgs84)

## Climate layers ##
chelsa = rast(chelsa_files)
# Crop first, then project onto the landscape grid
chelsa_cropped = terra::crop(chelsa, st_transform(extent_smaller, crs(spatial_data)))
chelsa_proj <- terra::project(chelsa_cropped, "EPSG:4326")
chelsa_proj = terra::project(chelsa_proj, spatial_data)
chelsa_proj <- mask(chelsa_proj, spatial_data[[1]])
spatial_data = c(spatial_data, chelsa_proj)

names(spatial_data) = c('elevation', 'pis', 'tree_canopy', 'ndvi', 'pd', 'wtd', 'hd', 'bio04', 'bio05', 'bio06', 'bio07', 'bio11', 'bio18', 'fcf', 'gdd10')

writeRaster(spatial_data, OUT_STACK, overwrite = T)

# Sites as SpatVector
sites_vect <- vect(site_coords, geom = c("lon", "lat"), crs = "EPSG:4326")

# Project to UTM zone 18N for buffering in metres
sites_proj <- project(sites_vect, "EPSG:32618")

# 500 m buffer around each site
buffers <- buffer(sites_proj, width = 500)

results <- site_coords

for (i in 1:nlyr(spatial_data)) {
  layer <- spatial_data[[i]]

  # Reproject raster to the buffer CRS
  if (!compareCRS(layer, buffers)) {
    layer <- project(layer, crs(buffers))
  }

  # Mean within buffer
  extracted <- terra::extract(layer, buffers, fun = mean, na.rm = TRUE)

  col_name <- names(spatial_data)[i]
  results[[col_name]] <- extracted[,2]
}

write.csv(results, OUT_TABLE, row.names = F)

# Values are written at most 11 characters wide (as in the original table,
# which was saved from a spreadsheet).
fmt11 <- function(x) vapply(x, function(v) {
  for (d in 10:1) { s <- format(signif(v, d), digits = d); if (nchar(s) <= 11) return(s) }
  s
}, character(1))
kept <- results[, c('site', 'lat', 'lon', KEPT)]
kept[KEPT] <- lapply(kept[KEPT], fmt11)
write.csv(kept, OUT_KEPT, row.names = F, quote = F)
say("Wrote %s, %s and %s", OUT_STACK, OUT_TABLE, OUT_KEPT)
