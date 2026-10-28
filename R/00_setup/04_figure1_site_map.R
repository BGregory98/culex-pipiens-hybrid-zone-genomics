# =============================================================================
# R/00_setup/04_figure1_site_map.R
# =============================================================================
# Maps the 24 field sites, numbered 1-24 as in Table S1, over the two
# environmental gradients: (A) percent impervious surface and (B) maximum
# summer temperature (bio05), with state and county boundaries, and writes the
# two panels side by side as Figure 1.
#
# Manuscript: Figure 1
#
# Inputs:  data/local_env_layers_2.tif                     (R/00_setup/01_site_predictors.R)
#          data/site_coords_with_buffered_spatial_data.csv (site coordinates; no script)
#          deposit: data/border_shapefiles/cb_2024_us_state_500k.shp
#          deposit: data/border_shapefiles/cb_2024_us_county_500k.shp
#                   US Census Bureau 2024 cartographic boundary files, 1:500,000
#                   (public; https://www.census.gov/geographies/mapping-files/
#                   time-series/geo/cartographic-boundary.html)
# Outputs: figures/final/Figure_1.png
#
# Runtime: under a minute.
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/sample_metadata.R")

library(terra)
library(tidyverse)
library(ggplot2)
library(tidyterra)
library(sf)
library(ggrepel)
library(patchwork)

RASTER_FILE <- 'data/local_env_layers_2.tif'
SITES_FILE  <- 'data/site_coords_with_buffered_spatial_data.csv'
need_files(RASTER_FILE, SITES_FILE)

OUT_FIG <- 'figures/final/Figure_1.png'
dir.create('figures/final', recursive = TRUE, showWarnings = FALSE)

states   <- st_read(deposit_file('data/border_shapefiles/cb_2024_us_state_500k.shp'),  quiet = TRUE)
counties <- st_read(deposit_file('data/border_shapefiles/cb_2024_us_county_500k.shp'), quiet = TRUE)

spatial_data <- rast(RASTER_FILE)

site_coords = read.csv(SITES_FILE)[,1:3]
site_coords = site_coords %>% arrange(site)
site_coords = site_coords %>% filter(site %in% field_names)

# sf object from lat/lon (coords = c(longitude, latitude)).
sites_sf <- st_as_sf(site_coords, coords = c("lon", "lat"), crs = 'EPSG:4326')

# Site numbers used on the map and in Table S1 (sites in alphabetical order).
sites_sf$site_num = c(24,  9,  8, 16,  2,  5, 17,  1, 10,  4, 18, 19,
                      21, 15,  3, 23, 11, 13,  7, 12, 22,  6, 20, 14)

# The figure is 13.33 in square at 300 dpi; text, points and the colour bar are
# sized for that canvas (0.85 x the sizes of a single 9 in panel).
site_map <- function(layer, fill_scale, legend_title, tag) {
  ggplot() +
    geom_spatraster(data = spatial_data[[layer]]) +
    fill_scale +
    geom_sf(data = states, fill = 'transparent', linewidth = 0.17) +
    geom_sf(data = counties, fill = 'transparent', linewidth = 0.17) +
    geom_sf(data = sites_sf, fill = 'deeppink', size = 5.1, shape = 21, color = 'black') +
    geom_text_repel(data = sites_sf, mapping = aes(label = site_num, geometry = geometry),
                    size = 5.1, stat = "sf_coordinates", point.size = 3.4, seed = 123) +
    annotate('text', x = -76.446, y = 39.490, label = tag,
             hjust = 1, vjust = 1, fontface = 'bold', size = 8.5) +
    coord_sf(xlim = c(-77.1824, -76.40722), ylim = c(38.69804, 39.51348), expand = F) +
    labs(fill = legend_title) +
    theme_bw(base_size = 15.3) +
    theme(axis.title            = element_blank(),
          legend.position       = 'bottom',
          legend.title.position = 'top',
          legend.title          = element_text(size = 18.7, hjust = 0.5),
          legend.text           = element_text(size = 20.4),
          legend.frame          = element_rect(colour = 'black'),
          legend.ticks          = element_blank()) +
    # Bar size set on the guide; a theme-level key width is stretched by patchwork.
    guides(fill = guide_colourbar(theme = theme(legend.key.width  = unit(3.27, 'in'),
                                                legend.key.height = unit(0.22, 'in'))))
}

# (A) Impervious surface
p_pis <- site_map(
  'pis',
  scale_fill_cross_blended_c(palette = 'warm_humid', limits = c(0, 90),
                             breaks = c(0, 45, 90), na.value = 'steelblue'),
  'Percent Impervious Surface', 'A')

# (B) bio05 (maximum summer temperature)
p_bio05 <- site_map(
  'bio05',
  scale_fill_distiller(palette = 'YlOrRd', direction = 1, limits = c(27.8, 32),
                       breaks = c(28, 30, 32), na.value = 'steelblue'),
  'Maximum Summer Temperature (°C)', 'B') +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig <- p_pis + p_bio05 + plot_layout(nrow = 1, guides = 'keep')

ggsave(OUT_FIG, fig, width = 4000 / 300, height = 4000 / 300, units = "in", dpi = 300)
say("Wrote %s", OUT_FIG)
