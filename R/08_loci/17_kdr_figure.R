# =============================================================================
# R/08_loci/17_kdr_figure.R
# =============================================================================
# Draws Figure S7, kdr (Vgsc L1014F) frequency across the 24 field pools:
#
#   A  frequency at each collection site (map, site numbers as Fig. 1)
#   B  frequency against maximum summer temperature      (n.s., P = 0.32)
#   C  frequency against percent impervious surface      (n.s., P = 0.16)
#   D  frequency against genome-wide PIP-like ancestry
#
# Panel D is included deliberately even though the correlation it shows does NOT
# survive a matched null (16_kdr_ancestry_matched_null.R, empirical P = 0.72).
# The caption says so. Showing the raw relationship and stating that it fails its
# null is more useful to a reader than omitting it, because the relationship is
# the one a reader would otherwise go looking for.
#
# One colour scale across all panels; the map background is neutral so the only
# hue in the figure encodes allele frequency.
#
# Manuscript: Figure S7
#
# Inputs:  results/resistance/Table_kdr_L1014F_by_site.csv    (14_kdr_frequency.R)
#          results/multilocus/per_pool_ancestry_per_pool.csv  (R/07_multi_locus/01_per_pool_ancestry.R)
#          deposit: data/border_shapefiles/cb_2024_us_state_500k.shp  (+ .shx/.dbf/.prj)
#          deposit: data/border_shapefiles/cb_2024_us_county_500k.shp (+ .shx/.dbf/.prj)
# Outputs: figures/final/supplement/Figure_S7.png
#
# Runtime: < 1 min
# =============================================================================

source("R/helpers/config.R")

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(sf); library(ggrepel); library(patchwork)
})

KDR_FILE <- 'results/resistance/Table_kdr_L1014F_by_site.csv'
ANC_FILE <- 'results/multilocus/per_pool_ancestry_per_pool.csv'
need_files(KDR_FILE, ANC_FILE)
STATE_SHP  <- deposit_file('data/border_shapefiles/cb_2024_us_state_500k.shp')
COUNTY_SHP <- deposit_file('data/border_shapefiles/cb_2024_us_county_500k.shp')

FIG_OUT <- 'figures/final/supplement/Figure_S7.png'
dir.create(dirname(FIG_OUT), recursive = TRUE, showWarnings = FALSE)

d <- fread(KDR_FILE)
anc <- fread(ANC_FILE)
d <- merge(d, anc[panel == 'ME1_private', .(site = pool, alpha_PIP = alpha_genome)], by = 'site')
stopifnot(nrow(d) == 24L)

SCALE <- scale_colour_viridis_c(name = 'kdr L1014F\nfrequency',
                                limits = c(0, 0.85), breaks = seq(0, 0.8, 0.2))
pts <- function(sz = 4) list(geom_point(aes(colour = kdr_L1014F), size = sz),
                             geom_point(shape = 21, fill = NA, colour = 'black',
                                        size = sz, stroke = 0.5))

# ---- A: map ------------------------------------------------------------------
states   <- st_read(STATE_SHP,  quiet = TRUE)
counties <- st_read(COUNTY_SHP, quiet = TRUE)
sf_sites <- st_as_sf(as.data.frame(d), coords = c('lon', 'lat'), crs = 'EPSG:4326')

pA <- ggplot() +
  geom_sf(data = states,   fill = 'grey99', colour = 'grey30', linewidth = 0.4) +
  geom_sf(data = counties, fill = NA,       colour = 'grey70', linewidth = 0.25) +
  geom_sf(data = sf_sites, aes(colour = kdr_L1014F), size = 4.5) +
  geom_sf(data = sf_sites, shape = 21, fill = NA, colour = 'black', size = 4.5, stroke = 0.5) +
  geom_text_repel(data = sf_sites, aes(label = site_num, geometry = geometry),
                  size = 3.4, stat = 'sf_coordinates', point.size = 6,
                  seed = 123, min.segment.length = 0.25) +
  SCALE +
  coord_sf(xlim = c(-77.1824, -76.40722), ylim = c(38.69804, 39.51348), expand = FALSE) +
  labs(tag = 'A') +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank(), panel.grid = element_line(colour = 'grey93'),
        axis.text = element_text(size = 8))

scat <- function(xv, xlab, tag) {
  ggplot(d, aes(.data[[xv]], kdr_L1014F)) +
    geom_smooth(method = 'lm', formula = y ~ x, colour = 'grey50',
                fill = 'grey88', linewidth = 0.6) +
    pts(3.6) + SCALE +
    labs(x = xlab, y = 'kdr L1014F frequency', tag = tag) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}
pB <- scat('bio05', 'maximum summer temperature (°C)', 'B')
pC <- scat('pis',   'percent impervious surface',           'C')
pD <- scat('alpha_PIP', 'genome-wide PIP-like ancestry', 'D')

fig <- (pA | (pB / pC / pD)) +
  plot_layout(widths = c(1.15, 1), guides = 'collect') &
  theme(legend.position = 'right')

ggsave(FIG_OUT, fig, width = 11, height = 9.5, units = 'in', dpi = 300)
cat(sprintf('wrote %s\n', FIG_OUT))
