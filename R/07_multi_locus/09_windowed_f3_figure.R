# =============================================================================
# R/07_multi_locus/09_windowed_f3_figure.R
# =============================================================================
# Figure S10: focal-window windowed f3 at CpL1a (top row, A-C) and CpL2a
# (bottom row, D-F) for the three reference-pair contrasts. For each field
# pool, the focal window's f3 z-score as a deviation from the median of that
# pool's genome-wide background (blue), with the background 2.5th-97.5th
# percentile range (grey whisker) centred at 0. Pools run left to right by
# maximum summer temperature (bio05, 500 m buffer mean; Table S1).
#
# Manuscript: Figure S10
#
# Inputs:  results/f3/focal_window_f3_cluster1.csv  (08_windowed_f3.R)
#          results/f3/focal_window_f3_cluster2.csv  (08_windowed_f3.R)
# Outputs: figures/final/supplement/Figure_S10.png
#
# Runtime: seconds
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/ref_colors.R")   # ref_pair_labels

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

IN_FILES <- sprintf("results/f3/focal_window_f3_%s.csv", c("cluster1", "cluster2"))
OUT_PNG  <- "figures/final/supplement/Figure_S10.png"
need_files(IN_FILES)

d <- rbindlist(lapply(IN_FILES, fread))

d[, locus := factor(c(cluster1 = "CpL1a", cluster2 = "CpL2a")[cluster], levels = c("CpL1a", "CpL2a"))]
d[, ref := factor(ref, levels = c("FL1-CAL1", "ME1-CAL1", "ME1-FL1"),
                  labels = unname(ref_pair_labels[c("FL1-CAL1", "ME1-CAL1", "ME1-FL1")]))]
pool_order <- unique(d[order(predictor), popA])
d[, popA := factor(popA, levels = pool_order)]
tags <- unique(d[, .(locus, ref)])[order(locus, ref)][, tag := LETTERS[.I]]

p <- ggplot(d, aes(popA)) +
  geom_errorbar(aes(ymin = bg_q025_z - bg_median_z, ymax = bg_q975_z - bg_median_z),
                colour = "grey60", width = 0.6) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_point(aes(y = 0), colour = "grey40", size = 1.2) +
  geom_point(aes(y = focal_z - bg_median_z), fill = "#377eb8", colour = "black", shape = 21, size = 2) +
  geom_text(data = tags, aes(x = -Inf, y = Inf, label = tag), inherit.aes = FALSE,
            hjust = -0.4, vjust = 1.3, fontface = "bold", size = 6) +
  facet_grid(locus ~ ref) +
  labs(x = "Field population", y = "f3 z-score (focal window vs.\ngenomic background)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_rect(fill = "grey85"),
        panel.grid.minor = element_blank())

dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)
ggsave(OUT_PNG, p, width = 10, height = 6.6, dpi = 300, bg = "white")
say("wrote %s", OUT_PNG)
