# =============================================================================
# R/07_multi_locus/05_linkage_decay_figure.R
# =============================================================================
# Draws Figure 4 from the binned LD summary written by 04_linkage_decay.R:
# median r^2 (interquartile band) against pairwise distance for SNP sets A-C
# and E, with the cross-locus unlinked set D as a dashed line, one facet per
# locus.
#
# Manuscript: Figure 4
#
# Inputs:  results/LD/Table_multi_locus_LD_decay.csv  (04_linkage_decay.R)
# Outputs: figures/final/Figure_4.png
#
# Runtime: seconds
# =============================================================================

source("R/helpers/config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

TABLE_IN <- 'results/LD/Table_multi_locus_LD_decay.csv'
FIG_OUT  <- 'figures/final/Figure_4.png'
need_files(TABLE_IN)

decay_summary <- read.csv(TABLE_IN, stringsAsFactors = FALSE)

# Locus display names for the locus keys in the table.
locus_labels <- c(chr1 = 'CpL1a', chr2 = 'CpL2a', chr3 = 'CpL3a')

# Bin midpoints for the distance bins of 04_linkage_decay.R.
bin_mid <- c('<1kb'=500,'1-5kb'=2.2e3,'5-20kb'=1e4,'20-100kb'=4.5e4,
             '100-250kb'=1.6e5,'250-500kb'=3.5e5,'500-750kb'=6.1e5,
             '750kb-1Mb'=8.7e5,'1-1.5Mb'=1.22e6,'1.5-2Mb'=1.73e6,'2-3Mb'=2.45e6,
             '3-5Mb'=3.9e6,'5-8Mb'=6.3e6,'8-12Mb'=9.8e6,'12-20Mb'=1.55e7,'>20Mb'=3e7)

MIN_PAIRS <- 10L

decay_plot <- decay_summary %>%
  filter(distance_bin != 'unlinked', n_pairs >= MIN_PAIRS) %>%
  mutate(distance_mid = bin_mid[as.character(distance_bin)],
         locus_label  = factor(locus_labels[locus], levels = unname(locus_labels)))

# Cross-locus unlinked appears as a horizontal dashed line (per-locus value)
unlinked_band <- decay_summary %>%
  filter(grepl('cross-locus unlinked', set)) %>%
  mutate(locus_label = factor(locus_labels[locus], levels = unname(locus_labels)))

set_colors <- c('A: outliers'             = '#AA4A44',
                'B: all-SNPs'             = '#4D72B0',
                'C: cluster-flank'        = '#E69F00',
                'D: cross-locus unlinked' = '#1f77b4',
                'E: chr baseline'         = 'grey40')

p_decay <- ggplot(decay_plot,
                  aes(x = distance_mid, y = median_r2,
                      color = set, fill = set)) +
  geom_hline(data = unlinked_band,
             aes(yintercept = median_r2, color = set),
             linetype = 'dashed', linewidth = 0.7) +
  geom_ribbon(aes(ymin = q25_r2, ymax = q75_r2),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  scale_x_log10(breaks = c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8),
                labels = c('1kb', '10kb', '100kb', '1Mb', '10Mb', '100Mb')) +
  scale_color_manual(values = set_colors, name = NULL) +
  scale_fill_manual(values  = set_colors, name = NULL) +
  # ggplot2 4.x renders the fill and colour scales as two separate legends here;
  # the published figure has one. Suppressing the fill guide merges them back.
  guides(fill = 'none') +
  facet_wrap(~ locus_label, ncol = 2, scales = 'fixed') +
  labs(x = 'Pairwise distance (log scale)',
       y = expression('Median ' * r^2 * ' across 24 field populations')) +
  theme_bw(base_size = 14) +
  theme(legend.position  = 'bottom',
        strip.background = element_blank(),
        strip.text       = element_text(face = 'bold', size = 12))

dir.create(dirname(FIG_OUT), showWarnings = FALSE, recursive = TRUE)
ggsave(FIG_OUT, p_decay, width = 10, height = 9, units = 'in', dpi = 300)
say('  Wrote %s', FIG_OUT)
