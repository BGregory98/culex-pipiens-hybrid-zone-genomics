# =============================================================================
# R/00_setup/06_temperature_vs_impervious.R
# =============================================================================
# Figure S1: maximum summer temperature (bio05) against percent impervious
# surface for the 24 field sites (Table S1 values), coloured by metropolitan
# area, with the Pearson correlation printed on the plot.
#
# Manuscript: Figure S1; Results "Mosquito collection and environmental
#             variation across the corridor" (r = 0.26); Supplementary Methods
#             "Environmental predictors"
#
# Inputs:  results/multilocus/Table_S1_predictors_full_precision.csv  (R/00_setup/05_table_s1_predictors.R)
# Outputs: figures/final/supplement/Figure_S1.png
#
# Runtime: seconds
# =============================================================================

source("R/helpers/config.R")
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

need_files("results/multilocus/Table_S1_predictors_full_precision.csv")

d <- fread("results/multilocus/Table_S1_predictors_full_precision.csv",
           showProgress = FALSE)
d <- d[grepl("^field", class)]

r_val <- cor(d$bio05, d$pis, use = "complete.obs")
say("r(bio05, impervious surface) = %.3f across %d sites", r_val, nrow(d))

metro_colors <- c(Baltimore = "#e66101", Washington = "#5e3c99")

p <- ggplot(d, aes(x = pis, y = bio05, colour = metro)) +
  geom_point(size = 2.5) +
  geom_text_repel(aes(label = site), size = 3.5,
                  max.overlaps = 30,
                  segment.colour = "grey60", segment.size = 0.3,
                  box.padding = 0.35, point.padding = 0.25,
                  seed = 42, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, label = sprintf("italic(r) == %.2f", r_val),
           parse = TRUE, hjust = -0.3, vjust = 1.6, size = 4.5) +
  scale_colour_manual(values = metro_colors, name = "Metropolitan area") +
  labs(x = "Percent impervious surface (%)",
       y = "Maximum summer temperature, bio05 (\u00B0C)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_line(linewidth = 0.25),
        legend.position = "bottom")

outfile <- "figures/final/supplement/Figure_S1.png"
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
ggsave(outfile, p, width = 7, height = 5, dpi = 300, bg = "white")
say("wrote %s", outfile)
