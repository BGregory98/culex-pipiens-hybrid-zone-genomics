# =============================================================================
# R/07_multi_locus/11_windowed_pool_pca_figure.R
# =============================================================================
# Figure S11 and Figure 3 from the 500 kb windows of 10_windowed_pool_pca.R.
# Figure S11:
#   A  PC1 % variance per window
#   B  Local-PCA MDS axis 2 (lostruct-style), which isolates CpL3a
#   C  Per-pool PC1 scores, oriented so + = molestus-like, coloured by
#      percent impervious surface
# Figure 3: the BayPass candidate SNPs (A) over PC1 variance per window (B), so
# the genome scan and the environment-free scan can be read against each other.
# Also prints the summary of the CpL3a windows quoted in Supplementary Results
# and the genome-wide ranking of the windows near the three loci.
#
# Manuscript: Figure 3; Figure S11; Results "A genome scan identifies three
#             candidate regions with distinct architectures" (the 25 most
#             structured windows) and the CpL3a section; Methods "Windowed PCA
#             and inversion breakpoints"; Supplementary Results "A genome-wide
#             windowed PCA recovers the inversion without prior coordinates"
#
# Inputs:  results/multilocus/windowed_pca/windowed_pca_500kb.csv         (10_windowed_pool_pca.R)
#          results/multilocus/windowed_pca/windowed_pca_scores_500kb.csv  (10_windowed_pool_pca.R)
#          results/baypass/v4/baypass_outliers_clustered.csv              (R/03_baypass/01_baypass_outliers.R)
# Outputs: figures/final/supplement/Figure_S11.png
#          figures/final/Figure_3.png
#
# Runtime: seconds
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")

library(tidyverse)
library(patchwork)

IN_DIR  <- "results/multilocus/windowed_pca"
RES_CSV    <- file.path(IN_DIR, "windowed_pca_500kb.csv")
SCORES_CSV <- file.path(IN_DIR, "windowed_pca_scores_500kb.csv")
OUT_PNG <- "figures/final/supplement/Figure_S11.png"
need_files(RES_CSV, SCORES_CSV)

chr_lab <- c(chr1 = "Chromosome 1", chr2 = "Chromosome 2", chr3 = "Chromosome 3")

res    <- read_csv(RES_CSV, show_col_types = FALSE) %>%
  mutate(chr = chr_lab[chr])
scores <- read_csv(SCORES_CSV, show_col_types = FALSE) %>%
  mutate(chr = chr_lab[chr])

# Candidate loci (this study). CpL3a coincides with inversion 3q of
# Cebrian-Camison et al. (2026).
loci <- tribble(
  ~chr,           ~start,            ~end,            ~label,
  "Chromosome 1", CPL1A_START / 1e6, CPL1A_END / 1e6, "CpL1a",
  "Chromosome 2", CPL2A_START / 1e6, CPL2A_END / 1e6, "CpL2a",
  "Chromosome 3", CPL3A_START / 1e6, CPL3A_END / 1e6, "CpL3a (3q)"
)

bands <- list(
  geom_rect(data = loci, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#1f77b4", alpha = 0.12)
)

base_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"))

facet <- facet_grid(. ~ chr, scales = "free_x", space = "free_x")

p95 <- quantile(res$pc1_pct, 0.95)
ymax_a <- max(res$pc1_pct) * 1.12   # label height

pA <- ggplot(res, aes(mid_mb, pc1_pct)) +
  bands +
  geom_hline(yintercept = p95, linetype = "dashed", colour = "#d62728", linewidth = 0.4) +
  geom_point(size = 0.7, colour = "grey20") +
  geom_text(data = loci, aes(x = (start + end) / 2, y = ymax_a, label = label),
            inherit.aes = FALSE, size = 3, fontface = "italic", colour = "#1f77b4") +
  facet +
  scale_y_continuous(limits = c(0, ymax_a * 1.08), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "PC1 variance (%)") +
  base_theme

pB <- ggplot(res, aes(mid_mb, mds2)) +
  bands +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 0.7, colour = "grey20") +
  facet +
  labs(x = NULL, y = "Local PCA\nMDS axis 2") +
  base_theme

pC <- ggplot(scores, aes(mid_mb, pc1_score, group = pool, colour = pis)) +
  bands +
  geom_line(linewidth = 0.3, alpha = 0.85) +
  scale_colour_viridis_c(name = "Impervious\nsurface (%)") +
  facet +
  labs(x = "Position (Mb)", y = "PC1 score\n(+ = molestus-like)") +
  base_theme +
  theme(legend.position = "bottom",
        legend.key.height = unit(0.35, "cm"),
        legend.key.width  = unit(1.2, "cm"))

fig <- pA / pB / pC + plot_annotation(tag_levels = "A")

dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)
ggsave(OUT_PNG, fig, width = 10, height = 7.5, dpi = 300, bg = "white")
say("wrote %s", OUT_PNG)


# ---- Summary of the CpL3a windows (500 kb) ----------------------------------

W <- 5e5
chr3 <- res %>% filter(chr == "Chromosome 3")
inside <- chr3 %>% filter(start >= CPL3A_START, start + W <= CPL3A_END)
near   <- chr3 %>% filter(start >= CPL3A_START - 2e6, start < CPL3A_END + 1.5e6)
far    <- res %>% filter(!(chr == "Chromosome 3" &
                           start >= CPL3A_START - 2e6 & start < CPL3A_END + 1.5e6))

say("\nWindowed PCA, 500 kb windows: %d windows", nrow(res))
say("Genome-wide median PC1 variance: %.1f%%", median(res$pc1_pct))
say("Windows entirely within CpL3a: %d (%.1f-%.1f Mb)",
    nrow(inside), min(inside$start) / 1e6, (max(inside$start) + W) / 1e6)
say("  PC1 variance: %.1f-%.1f%% (median %.1f%%)",
    min(inside$pc1_pct), max(inside$pc1_pct), median(inside$pc1_pct))
say("  MDS axis 2:   %.2f-%.2f", min(inside$mds2), max(inside$mds2))
say("  median r(PC1 score, impervious surface): %.2f", median(inside$cor_pis))
say("MDS axis 2 of chr3 windows from 2 Mb before to 1.5 Mb after CpL3a:")
for (i in seq_len(nrow(near))) {
  say("  %6.1f-%6.1f Mb  %5.2f", near$start[i] / 1e6, (near$start[i] + W) / 1e6, near$mds2[i])
}
say("Highest MDS axis 2 among all other windows: %.2f", max(far$mds2))


# ---- Ranking of the three loci in the environment-free scan (500 kb) ---------
# Which windows are most structured among the 24 field pools, with no
# environmental input: how far down the ranking the windows within 2 Mb of a
# candidate locus reach before the first window elsewhere.
BUF <- 2e6
loc_iv <- tribble(~chr, ~s, ~e,
  "Chromosome 1", CPL1A_START, CPL1A_END,
  "Chromosome 2", CPL2A_START, CPL2A_END,
  "Chromosome 3", CPL3A_START, CPL3A_END)
rk <- res %>%
  mutate(rank = rank(-pc1_pct, ties.method = "min"),
         near_locus = pmap_lgl(list(chr, start), function(ch, st) {
           iv <- loc_iv[loc_iv$chr == ch, ]
           st + W > iv$s - BUF & st < iv$e + BUF })) %>%
  arrange(rank)
first_out <- rk %>% filter(!near_locus) %>% slice(1)
top_in <- rk %>% filter(rank < first_out$rank)
say("\nRanking of all %d windows by PC1 variance (no environmental input):", nrow(rk))
say("  the %d most structured windows all lie within 2 Mb of a candidate locus",
    nrow(top_in))
say("    (PC1 %.1f-%.1f%%; CpL1a %d, CpL2a %d, CpL3a %d windows)",
    min(top_in$pc1_pct), max(top_in$pc1_pct),
    sum(top_in$chr == "Chromosome 1"), sum(top_in$chr == "Chromosome 2"),
    sum(top_in$chr == "Chromosome 3"))
say("  first window elsewhere: rank %d, %s %.2f Mb, PC1 %.1f%%",
    first_out$rank, first_out$chr, first_out$mid_mb, first_out$pc1_pct)
say("  genome-wide median %.1f%%, 95th percentile %.1f%%", median(res$pc1_pct), p95)

# the same ranking in 1 Mb windows
res1 <- read_csv(file.path(IN_DIR, "windowed_pca_1Mb.csv"), show_col_types = FALSE) %>%
  mutate(chr = chr_lab[chr])
rk1 <- res1 %>%
  mutate(rank = rank(-pc1_pct, ties.method = "min"),
         near_locus = pmap_lgl(list(chr, start), function(ch, st) {
           iv <- loc_iv[loc_iv$chr == ch, ]
           st + 1e6 > iv$s - BUF & st < iv$e + BUF })) %>%
  arrange(rank)
first1 <- rk1 %>% filter(!near_locus) %>% slice(1)
say("1 Mb windows: the %d most structured of %d lie within 2 Mb of a candidate locus; first elsewhere ranks %d (%s %.1f Mb, PC1 %.1f%%)",
    first1$rank - 1, nrow(rk1), first1$rank, first1$chr, first1$mid_mb, first1$pc1_pct)


# ---- Figure 3: the genome scan (A) over the environment-free scan (B) --------
# Panel A: the BayPass candidate SNPs (R/03_baypass/01_baypass_outliers.R),
# coloured by predictor category. Panel B: PC1 variance per 500 kb window from
# the windowed PCA above, which uses no environmental data. The three locus
# intervals are shaded in both, CpL2a at its 2.5 Mb block in each.
FOCAL_FILE <- "results/baypass/v4/baypass_outliers_clustered.csv"
OUT_FIG3   <- "figures/final/Figure_3.png"
need_files(FOCAL_FILE)
focal <- read_csv(FOCAL_FILE, show_col_types = FALSE)
chrom_levels <- sort(unique(focal$chromosome))
focal <- focal %>% mutate(chr = chr_lab[paste0("chr", match(chromosome, chrom_levels))])

bands3 <- tribble(
  ~locus,  ~chr,           ~start,           ~end,
  "CpL1a", "Chromosome 1", CPL1A_START,      CPL1A_END,
  "CpL2a", "Chromosome 2", CPL2A_START,      CPL2A_END,
  "CpL3a", "Chromosome 3", CPL3A_START,      CPL3A_END)
# common x extent per chromosome, so the two panels align
extent <- bind_rows(res %>% transmute(chr, x = 0),
                    res %>% transmute(chr, x = (start + W) / 1e6))

pA3 <- ggplot(focal, aes(x = position / 1e6, y = XtXst)) +
  geom_blank(data = extent, aes(x = x), inherit.aes = FALSE) +
  geom_rect(data = bands3, inherit.aes = FALSE,
            aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf),
            fill = "#1f77b4", alpha = 0.12) +
  geom_point(aes(colour = predictor_category), size = 0.7, alpha = 0.6) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = c(temperature = "brown1", urbanization = "cornflowerblue",
                                 other = "grey55"),
                      labels = c(temperature = "Temperature", urbanization = "Urbanization",
                                 other = "Other"),
                      breaks = c("temperature", "urbanization", "other"),
                      name = "Predictor category") +
  geom_text(data = bands3 %>% filter(locus != "CpL2a"), inherit.aes = FALSE,
            aes(x = start / 1e6, y = Inf, label = locus),
            hjust = 1, vjust = 2.2, nudge_x = -1, size = 3.5, fontface = "italic",
            colour = "#1f77b4") +
  geom_text(data = bands3 %>% filter(locus == "CpL2a"), inherit.aes = FALSE,
            aes(x = end / 1e6, y = Inf, label = locus),
            hjust = 0, vjust = 2.2, nudge_x = 1, size = 3.5, fontface = "italic",
            colour = "#1f77b4") +
  labs(x = NULL, y = expression("XtX"[st])) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", legend.margin = margin(t = -5),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 13))

pB3 <- ggplot(res, aes(mid_mb, pc1_pct)) +
  geom_blank(data = extent, aes(x = x), inherit.aes = FALSE) +
  geom_rect(data = loci, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#1f77b4", alpha = 0.12) +
  geom_hline(yintercept = p95, linetype = "dashed", colour = "#d62728", linewidth = 0.4) +
  geom_point(size = 0.8, colour = "grey20") +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  labs(x = "Position (Mb)", y = "PC1 variance (%)") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(), strip.text = element_blank())

fig3 <- pA3 / pB3 + plot_layout(heights = c(1.6, 1)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold"))
ggsave(OUT_FIG3, fig3, width = 12, height = 7.5, units = "in", dpi = 300, bg = "white")
say("wrote %s", OUT_FIG3)
