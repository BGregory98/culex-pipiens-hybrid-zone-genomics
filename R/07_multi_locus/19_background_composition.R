# =============================================================================
# R/07_multi_locus/19_background_composition.R
# =============================================================================
# Three analyses that share one pass over the allele-frequency table.
#   1. The differentiation-matched null of 14_matched_null_xtx.R, further
#      restricted to regions whose panel holds between half and twice the focal
#      panel's SNPs (small panels give shifted correlations).
#   2. Background composition from private-panel alpha outside all three
#      candidate intervals, at three panel thresholds. FST is deflated by
#      within-population diversity and the three references rank-order
#      inversely by their own heterozygosity, so the background is re-estimated
#      from the diversity-independent private panels.
#   3. Stability of the six-SNP QUI-private panel at CpL2a: Delta slope at three
#      thresholds, the six per-SNP series, leave-one-out, and a bootstrap over
#      the six SNPs.
#
# Manuscript: Supplementary Results "Locus statistics against a differentiation-
#             matched null" (panel-size-matched regions), "Background composition
#             estimated from lineage-private panels" and "Stability of the
#             six-SNP QUI-private panel at CpL2a"; Results "Corridor populations
#             are weakly differentiated and predominantly Cx. p. f. pipiens-like"
#             (ordering holds across panel thresholds)
#
# Inputs:  results/multilocus/Table_xtx_matched_null_windows.csv   (R/07_multi_locus/14_matched_null_xtx.R)
#          results/multilocus/Table_xtx_matched_null_summary.csv   (R/07_multi_locus/14_matched_null_xtx.R)
#          results/multilocus/Table_xtx_matched_null_CpL2a_QUI_windows.csv, _summary.csv
#                                                                  (R/07_multi_locus/15_matched_null_cpl2a_qui.R)
#          results/multilocus/Table_per_pool_diversity.csv         (R/04_pop_structure/05_pool_diversity.R)
#          data/genome_wide_FST_for_field_vs_references.csv        (R/04_pop_structure/01_pairwise_fst.R)
#          data/site_coords_with_spatial_data_allvars.csv          (R/00_setup/01_site_predictors.R)
#          deposit: data/allele_frequencies_v4.csv                 (R/01_filtering/01_allele_frequency_table.R)
# Outputs: results/multilocus/Table_background_composition.csv
#
# Runtime: under a minute
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")
suppressMessages(library(data.table))

AF_PATH <- deposit_file("data/allele_frequencies_v4.csv")
need_files("results/multilocus/Table_xtx_matched_null_windows.csv",
           "results/multilocus/Table_xtx_matched_null_summary.csv",
           "results/multilocus/Table_xtx_matched_null_CpL2a_QUI_windows.csv",
           "results/multilocus/Table_xtx_matched_null_CpL2a_QUI_summary.csv",
           "results/multilocus/Table_per_pool_diversity.csv",
           "data/genome_wide_FST_for_field_vs_references.csv",
           "data/site_coords_with_spatial_data_allvars.csv")

set.seed(42)   # the bootstrap over SNPs in part 3

THRESH <- list(c(0.95, 0.05), c(0.90, 0.10), c(0.85, 0.15))
# CpL2a is the 500 kb candidate core (13.5-14.0 Mb) in the private-panel analyses.
LOCI <- list(
  CpL1a = list(chr = CHR1, start = CPL1A_START, end = CPL1A_END, env = "bio05"),
  CpL2a = list(chr = CHR2, start = CPL2A_CORE_START,      end = CPL2A_CORE_END,    env = "bio05"),
  CpL3a = list(chr = CHR3, start = CPL3A_START, end = CPL3A_END, env = "pis")
)

# -----------------------------------------------------------------------------
# 1. Differentiation-matched null, additionally matched on panel size
# -----------------------------------------------------------------------------
# The three donor panels from 14, plus the CpL2a QUI-private panel from 15.
w <- rbind(fread("results/multilocus/Table_xtx_matched_null_windows.csv", showProgress = FALSE),
           fread("results/multilocus/Table_xtx_matched_null_CpL2a_QUI_windows.csv",
                 showProgress = FALSE)[locus == "CpL2a_QUI"])
s <- rbind(fread("results/multilocus/Table_xtx_matched_null_summary.csv", showProgress = FALSE),
           fread("results/multilocus/Table_xtx_matched_null_CpL2a_QUI_summary.csv",
                 showProgress = FALSE)[locus == "CpL2a_QUI"])
say("=== differentiation-matched null, additionally matched on panel size ===")
say("P(dir) is the share of matched regions at least as extreme as the focal r in its own direction.")
say("%-9s %8s %7s %9s %9s %9s %11s %11s", "locus", "focal_n", "matched", "null mean",
    "null med", "focal r", "P(|r|>=)", "P(dir)")
for (L in unique(w$locus)) {
  d  <- w[locus == L & is.finite(r)]
  fo <- s[locus == L]
  k  <- d[n_panel >= fo$focal_n_panel / 2 & n_panel <= fo$focal_n_panel * 2]
  if (nrow(k) < 20) {
    say("%-9s %8d %7d   too few panel-size-matched regions to form a null",
        L, fo$focal_n_panel, nrow(k))
    say("        (focal panel is enriched relative to comparable regions -- reported as such)")
    next
  }
  p_dir <- if (fo$focal_r >= 0) mean(k$r >= fo$focal_r) else mean(k$r <= fo$focal_r)
  say("%-9s %8d %7d %9.3f %9.3f %9.3f %11.4f %11.4f", L, fo$focal_n_panel, nrow(k),
      mean(k$r), median(k$r), fo$focal_r,
      mean(abs(k$r) >= abs(fo$focal_r)), p_dir)
}

# -----------------------------------------------------------------------------
# 2. Read the matrix once; genome-wide panels outside all candidate intervals
# -----------------------------------------------------------------------------
hdr   <- names(fread(AF_PATH, nrows = 0, showProgress = FALSE))
pools <- setdiff(hdr, c("chromosome", "position", "cumulative_position"))
field <- setdiff(pools, REFS)
env   <- fread("data/site_coords_with_spatial_data_allvars.csv", showProgress = FALSE)
setkey(env, site); e <- env[field]

say("")
say("reading %s ...", AF_PATH)
af <- fread(AF_PATH, select = c("chromosome", "position", pools), showProgress = FALSE)
say("matrix rows: %s", format(nrow(af), big.mark = ","))

q <- af[[QUI]]; p <- af[[PIP]]; m <- af[[MOL]]
ok <- is.finite(q) & is.finite(p) & is.finite(m)
in_cand <- rep(FALSE, nrow(af))
for (L in LOCI) in_cand <- in_cand | (af$chromosome == L$chr &
                                      af$position >= L$start & af$position <= L$end)
F <- as.matrix(af[, ..field])

alpha_of <- function(mask) if (sum(mask) < 3) rep(NA_real_, length(field)) else
  colMeans(F[mask, , drop = FALSE], na.rm = TRUE)

say("")
say("=== background composition from private-panel alpha, outside all candidate intervals ===")
bg <- list()
for (th in THRESH) {
  HI <- th[1]; LO <- th[2]
  masks <- list(
    QUI = ok & !in_cand & q >= HI & p <= LO & m <= LO,
    PIP = ok & !in_cand & p >= HI & q <= LO & m <= LO,
    MOL = ok & !in_cand & m >= HI & q <= LO & p <= LO)
  a <- lapply(masks, alpha_of)
  say("  thresholds >= %.2f / <= %.2f : n = %s QUI, %s PIP, %s MOL", HI, LO,
      format(sum(masks$QUI), big.mark = ","), format(sum(masks$PIP), big.mark = ","),
      format(sum(masks$MOL), big.mark = ","))
  say("    mean alpha across 24 field pools: QUI %.3f  PIP %.3f  MOL %.3f",
      mean(a$QUI), mean(a$PIP), mean(a$MOL))
  say("    range QUI %.3f-%.3f | PIP %.3f-%.3f | MOL %.3f-%.3f",
      min(a$QUI), max(a$QUI), min(a$PIP), max(a$PIP), min(a$MOL), max(a$MOL))
  if (identical(th, THRESH[[1]])) bg <- a
}
say("  NOTE: alpha values are indices, not proportions, and are not expected to sum to 1.")

# The same at >= 0.95 / <= 0.05 with the candidate intervals included: the
# genome-wide panels of Table S4 (Results, background composition).
gw_masks <- list(QUI = ok & q >= 0.95 & p <= 0.05 & m <= 0.05,
                 PIP = ok & p >= 0.95 & q <= 0.05 & m <= 0.05,
                 MOL = ok & m >= 0.95 & q <= 0.05 & p <= 0.05)
gw_a <- lapply(gw_masks, alpha_of)
say("  whole genome (intervals included): n = %s QUI, %s PIP, %s MOL; mean alpha QUI %.3f  PIP %.3f  MOL %.3f",
    format(sum(gw_masks$QUI), big.mark = ","), format(sum(gw_masks$PIP), big.mark = ","),
    format(sum(gw_masks$MOL), big.mark = ","), mean(gw_a$QUI), mean(gw_a$PIP), mean(gw_a$MOL))

say("")
say("=== is the FST ordering explained by reference heterozygosity? ===")
het <- fread("results/multilocus/Table_per_pool_diversity.csv", showProgress = FALSE)
gw  <- fread("data/genome_wide_FST_for_field_vs_references.csv", showProgress = FALSE)
fs  <- gw[, .(meanFST = mean(gw_fst)), by = ref]
fs[, lineage := c(CAL1 = "MOL", FL1 = "QUI", ME1 = "PIP")[ref]]
fs[, refhet := het$exp_het[match(ref, het$pool)]]
fs[, bgalpha := sapply(lineage, function(l) mean(bg[[l]]))]
print(fs[order(meanFST), .(lineage, meanFST = round(meanFST, 3),
                           ref_het = round(refhet, 3), bg_alpha = round(bgalpha, 3))],
      row.names = FALSE)
say("  Spearman(meanFST, reference heterozygosity) = %+.2f",
    cor(fs$meanFST, fs$refhet, method = "spearman"))
say("  -> FST rank-orders inversely with the reference's own diversity, so the")
say("     'least differentiated from PIP' ordering is not by itself an ancestry result.")

# -----------------------------------------------------------------------------
# 3. CpL2a QUI-private panel: Delta sensitivity, per-SNP series, bootstrap
# -----------------------------------------------------------------------------
L <- LOCI$CpL2a
iv <- af$chromosome == L$chr & af$position >= L$start & af$position <= L$end
x  <- as.numeric(e[[L$env]])
say("")
say("=== CpL2a QUI-private panel: threshold sensitivity on Delta ===")
say("%-16s %6s %10s %9s %11s", "thresholds", "n", "slope/C", "adjR2", "P")
for (th in THRESH) {
  HI <- th[1]; LO <- th[2]
  loc <- ok & iv & q >= HI & p <= LO & m <= LO
  gen <- ok & !in_cand & q >= HI & p <= LO & m <= LO
  if (sum(loc) < 3) { say("%-16s %6d   too few", sprintf(">=%.2f/<=%.2f", HI, LO), sum(loc)); next }
  del <- alpha_of(loc) - alpha_of(gen)
  f <- lm(del ~ x); ss <- summary(f)
  say("%-16s %6d %+10.4f %9.3f %11.3g", sprintf(">=%.2f/<=%.2f", HI, LO),
      sum(loc), coef(f)[2], ss$adj.r.squared, ss$coefficients[2, 4])
}

say("")
say("=== the six QUI-private SNPs at CpL2a, one series each ===")
sel <- which(ok & iv & q >= 0.95 & p <= 0.05 & m <= 0.05)
say("%-14s %8s %8s %10s %11s", "position", "cool", "warm", "slope/C", "P")
for (i in sel) {
  y <- as.numeric(F[i, ])
  f <- lm(y ~ x); ss <- summary(f)
  say("%-14s %8.3f %8.3f %+10.4f %11.3g", format(af$position[i], big.mark = ","),
      mean(y[order(x)[1:6]]), mean(y[order(x, decreasing = TRUE)[1:6]]),
      coef(f)[2], ss$coefficients[2, 4])
}
Y <- F[sel, , drop = FALSE]
say("  leave-one-out slope of the panel mean (drop each SNP in turn):")
for (j in seq_len(nrow(Y))) {
  a <- colMeans(Y[-j, , drop = FALSE])
  say("    without SNP %d: slope %+.4f  P = %.3g", j,
      coef(lm(a ~ x))[2], summary(lm(a ~ x))$coefficients[2, 4])
}

say("")
say("=== bootstrap over the six SNPs (1,000 resamples) ===")
sl <- replicate(1000, {
  k <- sample(nrow(Y), nrow(Y), replace = TRUE)
  coef(lm(colMeans(Y[k, , drop = FALSE]) ~ x))[2]
})
say("  observed slope %+.4f ; bootstrap 95%% interval %+.4f to %+.4f ; median %+.4f",
    coef(lm(colMeans(Y) ~ x))[2], quantile(sl, .025), quantile(sl, .975), median(sl))
say("  resamples with a negative slope: %.1f%%", 100 * mean(sl < 0))

fwrite(data.table(lineage = fs$lineage, mean_fst = fs$meanFST,
                  ref_het = fs$refhet, bg_alpha = fs$bgalpha),
       "results/multilocus/Table_background_composition.csv")
say("")
say("wrote results/multilocus/Table_background_composition.csv")
