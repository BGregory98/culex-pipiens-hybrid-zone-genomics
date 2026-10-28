# =============================================================================
# R/07_multi_locus/23_locus_vs_background_ancestry.R
# =============================================================================
# Whether each locus's environmental association survives when the locus is
# allowed to follow genome-wide donor ancestry at its own rate.
#
# Delta = alpha_locus - alpha_genome assumes that a locus gains donor ancestry
# one-for-one with the genome. That holds when the donor-private panel marks
# donor ancestry equally well inside the locus and across the genome. It need
# not hold at an inversion: the panels are defined from single colonies, so
# genome-wide many panel alleles are imperfect markers of wild donor ancestry,
# whereas suppressed recombination between arrangements preserves fixed
# differences inside an inversion. A better-marked locus then rises faster than
# the background as donor ancestry rises, and Delta inherits any environmental
# trend in the background with no locus-specific effect.
#
# Per locus: alpha_locus ~ alpha_genome (the slope, 1 under the Delta
# assumption), and the joint model alpha_locus ~ alpha_genome + predictor on the
# linear and log scales. The joint model is the comparison that matters: it asks
# which of the two, genome-wide ancestry or the environment, predicts the locus
# once the other is known. The slope on alpha_genome alone is not a measure of
# marker quality where the predictor drives both the locus and the background
# (at CpL1a it is 12.6), so it is read only alongside the joint model.
#
# Manuscript: Results "A molestus inversion carries urban molestus-like
#             ancestry, independently of temperature"; Methods "Per-population
#             donor-ancestry test"; Supplementary Methods "Donor-private ancestry
#             panels and the per-population donor-ancestry test"; Supplementary
#             Results "The CpL3a urban association against genome-wide
#             molestus-like ancestry"; Discussion (CpL3a)
#
# Inputs:  results/multilocus/per_pool_ancestry_per_pool.csv  (R/07_multi_locus/01_per_pool_ancestry.R)
# Outputs: results/multilocus/Table_locus_vs_background_ancestry.csv
#
# Runtime: seconds
# =============================================================================

source("R/helpers/config.R")
suppressMessages(library(data.table))

need_files("results/multilocus/per_pool_ancestry_per_pool.csv")
pp <- fread("results/multilocus/per_pool_ancestry_per_pool.csv")

# The three tested panels (the 3-SNP chr2 ME1-private panel is not tested).
TESTS <- data.table(
  locus = c("CpL1a", "CpL2a", "CpL3a"),
  lkey  = c("chr1", "chr2_nonFL1_deficit", "chr3"),
  panel = c("FL1_private", "non_FL1_deficit", "CAL1_private"))

out <- rbindlist(lapply(seq_len(nrow(TESTS)), function(i) {
  t <- TESTS[i]
  x <- pp[locus == t$lkey & panel == t$panel]
  stopifnot(nrow(x) == 24)
  f0 <- lm(alpha_focal ~ alpha_genome, data = x)
  f1 <- lm(alpha_focal ~ alpha_genome + env, data = x)
  fl <- if (all(x$alpha_focal > 0 & x$alpha_genome > 0))
          lm(log(alpha_focal) ~ log(alpha_genome) + env, data = x) else NULL
  fd <- lm(delta ~ env, data = x)
  s0 <- summary(f0); s1 <- summary(f1)
  data.table(locus = t$locus, panel = t$panel, n = nrow(x),
             slope_on_background = coef(f0)[2], slope_se = coef(s0)[2, 2],
             r2_on_background = s0$r.squared,
             delta_env_p_onesided = pt(coef(summary(fd))[2, 3], df = fd$df.residual,
                                       lower.tail = FALSE),
             env_p_given_background = coef(s1)["env", 4],
             background_p_given_env = coef(s1)["alpha_genome", 4],
             env_p_given_background_log = if (is.null(fl)) NA_real_ else coef(summary(fl))["env", 4],
             log_slope_on_background = if (is.null(fl)) NA_real_ else coef(fl)[2],
             r_background_env = cor(x$alpha_genome, x$env))
}))
fwrite(out, "results/multilocus/Table_locus_vs_background_ancestry.csv")

say("Locus ancestry against genome-wide ancestry (24 field pools)")
for (i in seq_len(nrow(out))) {
  o <- out[i]
  say("")
  say("%s (%s)", o$locus, o$panel)
  say("  alpha_locus ~ alpha_genome: slope %.2f (SE %.2f), R2 %.2f   [Delta assumes 1]",
      o$slope_on_background, o$slope_se, o$r2_on_background)
  say("  predictor given background: P = %.3g (linear), %.3g (log scale; log-log slope %.2f)",
      o$env_p_given_background, o$env_p_given_background_log, o$log_slope_on_background)
  say("  background given predictor: P = %.3g", o$background_p_given_env)
  say("  Delta ~ predictor, one-sided P = %.3g", o$delta_env_p_onesided)
  say("  r(genome-wide alpha, predictor) = %+.2f", o$r_background_env)
}
say("\nwrote results/multilocus/Table_locus_vs_background_ancestry.csv")
