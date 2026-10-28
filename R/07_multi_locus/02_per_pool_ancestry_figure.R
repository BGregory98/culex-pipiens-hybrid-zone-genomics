# =============================================================================
# R/07_multi_locus/02_per_pool_ancestry_figure.R
# =============================================================================
# Draws Figure 5 from the tables written by 01_per_pool_ancestry.R: per-pool
# Delta against the locus's predictor at CpL1a (QUI panel, bio05), CpL2a
# (non-QUI deficit panel, bio05) and CpL3a (MOL panel, percent impervious
# surface), with bootstrap intervals, an OLS fit, and the Wilcoxon and
# slope-permutation P values as panel subtitles.
#
# Manuscript: Figure 5
#
# Inputs:  results/multilocus/per_pool_ancestry_per_pool.csv  (01_per_pool_ancestry.R)
#          results/multilocus/per_pool_ancestry_tests.csv     (01_per_pool_ancestry.R)
# Outputs: figures/final/Figure_5.png
#
# Runtime: seconds
# =============================================================================

source("R/helpers/config.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

PER_POOL <- 'results/multilocus/per_pool_ancestry_per_pool.csv'
TESTS    <- 'results/multilocus/per_pool_ancestry_tests.csv'
FIG_OUT  <- 'figures/final/Figure_5.png'
need_files(PER_POOL, TESTS)

per_pool <- read.csv(PER_POOL, stringsAsFactors = FALSE)
tests    <- read.csv(TESTS,    stringsAsFactors = FALSE)

# ---- panel definitions ------------------------------------------------------
main_panels <- c(chr1 = 'FL1_private',
                 chr2 = 'non_FL1_deficit',
                 chr3 = 'CAL1_private')

# The deficit panel's rows were written under locus_label 'chr2_nonFL1_deficit',
# not 'chr2'. Cosmetic lookups stay keyed by chromosome; the data filter uses
# the stored locus label.
main_loci <- c(chr1 = 'chr1',
               chr2 = 'chr2_nonFL1_deficit',
               chr3 = 'chr3')

# Panel titles use the lineage display labels (FL1 = QUI, CAL1 = MOL).
panel_titles <- c(
  chr1 = expression('CpL1a'~~italic(alpha)['QUI']),
  chr2 = expression('CpL2a'~~italic(alpha)['non-QUI']),
  chr3 = expression('CpL3a'~~italic(alpha)['MOL'])
)

panel_env_labels <- c(
  chr1 = expression('bio05 (max temp.,'~degree*C*')'),
  chr2 = expression('bio05 (max temp.,'~degree*C*')'),
  chr3 = '% impervious surface'
)

make_panel <- function(locus_key) {
  df <- per_pool[per_pool$locus == main_loci[locus_key] &
                 per_pool$panel == main_panels[locus_key], ]
  test_row <- tests[tests$locus == main_loci[locus_key] &
                    tests$panel == main_panels[locus_key], ]

  subtitle <- sprintf('Wilcoxon p = %s    perm p (slope) = %s',
                      formatC(test_row$wilcoxon_p,   format = 'g', digits = 2),
                      formatC(test_row$perm_p_slope, format = 'g', digits = 2))

  ggplot(df, aes(x = env, y = delta)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey50') +
    geom_smooth(method = 'lm', formula = y ~ x, color = 'firebrick',
                fill = 'firebrick', alpha = 0.15) +
    geom_errorbar(aes(ymin = delta_ci_lo, ymax = delta_ci_hi),
                  width = 0, color = 'grey40', alpha = 0.6) +
    geom_point(size = 2.5, color = 'black') +
    scale_y_continuous(limits = c(-0.3, 0.3)) +
    labs(title = panel_titles[[locus_key]],
         subtitle = subtitle,
         x = panel_env_labels[[locus_key]],
         y = expression(Delta[i]~~'='~~alpha[locus]-alpha[genome])) +
    theme_bw(base_size = 18) +
    theme(plot.title    = element_text(face = 'bold'),
          plot.subtitle = element_text(size = 9, color = 'grey30'))
}

fig <- make_panel('chr1') + make_panel('chr2') + make_panel('chr3') +
       plot_layout(ncol = 3)

dir.create(dirname(FIG_OUT), showWarnings = FALSE, recursive = TRUE)
ggsave(FIG_OUT, fig, width = 14, height = 6, dpi = 300)
say('  Wrote %s', FIG_OUT)
