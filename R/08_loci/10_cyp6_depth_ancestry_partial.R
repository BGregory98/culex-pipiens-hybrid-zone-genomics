# =============================================================================
# R/08_loci/10_cyp6_depth_ancestry_partial.R
# =============================================================================
# Two checks on the CpL1a cluster depth ratio, from cached tables:
#
#  1. Ancestry budget and reference ordering. The genome-wide private-panel
#     alphas sum to well below 1, whereas at CpL1a a "residual-as-PIP"
#     subtraction forces the three to sum to 1; and the reference cluster
#     ratios are compared with the ordering divergence-driven mismapping
#     predicts (PIP < MOL < QUI).
#  2. Partial correlation. Depth ~ alpha_QUI is offered as what places the
#     extra copies on the incoming material, but alpha_QUI itself tracks bio05,
#     so the two may be one relationship re-expressed. Reports
#     r(depth, alpha_QUI | bio05) and r(depth, bio05 | alpha_QUI).
#
# The calibration of the field depth cline against the three references is in
# 09_cyp6_reference_pool_depth.R (section 3), on the frequency of the CpL1a
# QUI-state haplotype. An earlier mixture prediction here weighted the colony
# depths by panel alpha, which is an index rather than a proportion, and treated
# PIP as a single copy state; it was removed (2026-10-03).
#
# Manuscript: Results "Cx. quinquefasciatus ancestry and sequencing read depth
#             at a CYP6 cluster rise with summer temperature" (partial
#             correlation +0.60); Supplementary Results "Whether the CpL1a depth
#             cline could be a mismapping artifact" (partial correlations,
#             reference ordering)
#
# Inputs:  data/site_coords_with_spatial_data_allvars.csv          (tracked input; no producing script in R/)
#          results/multilocus/Table_private_panel_alpha_per_pool.csv   (R/07_multi_locus/12_private_panel_clines.R)
#          results/multilocus/Table_background_composition.csv    (R/07_multi_locus/19_background_composition.R)
#          results/chr1/Table_cyp6_cluster_depth_with_references.csv   (09_cyp6_reference_pool_depth.R)
# Outputs: results/chr1/Table_cyp6_depth_partial.csv
#
# Runtime: < 1 min
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")

suppressMessages({
  library(data.table)
  if (requireNamespace("bit64", quietly = TRUE)) library(bit64)
})

ENV_FILE   <- "data/site_coords_with_spatial_data_allvars.csv"
ALPHA_FILE <- "results/multilocus/Table_private_panel_alpha_per_pool.csv"
BG_FILE    <- "results/multilocus/Table_background_composition.csv"
DEPTH_FILE <- "results/chr1/Table_cyp6_cluster_depth_with_references.csv"
need_files(ENV_FILE, ALPHA_FILE, BG_FILE, DEPTH_FILE)

env  <- fread(ENV_FILE, showProgress = FALSE)
al   <- fread(ALPHA_FILE, showProgress = FALSE)
bg   <- fread(BG_FILE, showProgress = FALSE)
dep  <- fread(DEPTH_FILE, showProgress = FALSE)

field <- field_names
setkey(env, site)
x_bio <- as.numeric(env[field][["bio05"]])

aq <- al[locus == "CpL1a" & panel == "QUI-private"]; setkey(aq, pool)
am <- al[locus == "CpL1a" & panel == "MOL-private"]; setkey(am, pool)
a_qui <- aq[field][["alpha"]]
a_mol <- am[field][["alpha"]]
setkey(dep, pool)
d_clu <- dep[field][["cluster_ratio"]]
stopifnot(!anyNA(a_qui), !anyNA(a_mol), !anyNA(d_clu))

r_QUI <- dep[pool == QUI, cluster_ratio]
r_PIP <- dep[pool == PIP, cluster_ratio]
r_MOL <- dep[pool == MOL, cluster_ratio]

# -----------------------------------------------------------------------------
# 1. Ancestry budget and reference ordering
# -----------------------------------------------------------------------------
say("=== 1. ancestry budget and reference ordering ===")
s <- sum(bg$bg_alpha)
say("genome-wide alpha sum: PIP %.3f + MOL %.3f + QUI %.3f = %.3f",
    bg[lineage == "PIP", bg_alpha], bg[lineage == "MOL", bg_alpha],
    bg[lineage == "QUI", bg_alpha], s)
say("  -> the three panels recover %.1f%% of a unit ancestry budget; %.1f%% unaccounted",
    100 * s, 100 * (1 - s))

o <- order(x_bio)
cool <- head(o, 6); warm <- tail(o, 6)
aq_c <- mean(a_qui[cool]); aq_w <- mean(a_qui[warm])
am_c <- mean(a_mol[cool]); am_w <- mean(a_mol[warm])
say("")
say("CpL1a panel means, coolest six vs warmest six:")
say("  alpha_QUI %.3f -> %.3f  (change %+.3f)", aq_c, aq_w, aq_w - aq_c)
say("  alpha_MOL %.3f -> %.3f  (change %+.3f)", am_c, am_w, am_w - am_c)
say("  residual-as-PIP %.3f -> %.3f   <-- the invalid subtraction", 1 - aq_c - am_c,
    1 - aq_w - am_w)
say("  at CpL1a the three 'sum' to 1.000 BY CONSTRUCTION, against %.3f genome-wide.", s)

say("")
say("reference cluster ratios: MOL %.3f, PIP %.3f, QUI %.3f", r_MOL, r_PIP, r_QUI)
say("  PIP-to-QUI gap = %.3f", r_QUI - r_PIP)
say("divergence from the Swedish PIP-lineage assembly runs PIP < MOL < QUI,")
say("so divergence-driven mismapping predicts that same depth ordering.")
say("  ordering predicted by divergence-driven mismapping : PIP < MOL < QUI")
say("  ordering observed  : MOL %.3f < PIP %.3f < QUI %.3f", r_MOL, r_PIP, r_QUI)
say("  -> MOL is the more divergent of the two Cx. pipiens references yet sits")
say("     %.3f BELOW PIP. The prediction fails at the discriminating comparison.",
    r_PIP - r_MOL)

# -----------------------------------------------------------------------------
# 2. Partial correlation
# -----------------------------------------------------------------------------
pcor <- function(a, b, ctrl) {
  ra <- residuals(lm(a ~ ctrl)); rb <- residuals(lm(b ~ ctrl))
  n <- length(a); r <- cor(ra, rb)
  tstat <- r * sqrt((n - 3) / (1 - r^2))
  list(r = r, p = 2 * pt(-abs(tstat), df = n - 3))
}
say("")
say("=== 2. is the depth-ancestry association independent of temperature? ===")
say("simple correlations:")
say("  r(depth, alpha_QUI) = %+.3f   (adj R2 %.3f)", cor(d_clu, a_qui),
    summary(lm(d_clu ~ a_qui))$adj.r.squared)
say("  r(depth, bio05)     = %+.3f   (adj R2 %.3f)", cor(d_clu, x_bio),
    summary(lm(d_clu ~ x_bio))$adj.r.squared)
say("  r(alpha_QUI, bio05) = %+.3f   (adj R2 %.3f)", cor(a_qui, x_bio),
    summary(lm(a_qui ~ x_bio))$adj.r.squared)
pa <- pcor(d_clu, a_qui, x_bio)
pb <- pcor(d_clu, x_bio, a_qui)
say("")
say("partial correlations:")
say("  r(depth, alpha_QUI | bio05) = %+.3f , P = %.4g",
    pa$r, pa$p)
say("  r(depth, bio05 | alpha_QUI) = %+.3f , P = %.4g", pb$r, pb$p)
f3 <- lm(d_clu ~ x_bio + a_qui)
say("  joint model: P(bio05) = %.4g , P(alpha_QUI) = %.4g",
    summary(f3)$coefficients[2, 4], summary(f3)$coefficients[3, 4])
say("")
if (pa$p < 0.05 && pb$p >= 0.05) {
  say("  -> the ancestry association SURVIVES the temperature control while the")
  say("     temperature association does not survive the ancestry control.")
} else if (pa$p < 0.05) {
  say("  -> both partials survive; the two predictors carry separable information.")
} else {
  say("  -> the ancestry partial does NOT survive. Report that both track")
  say("     temperature and that haplotype linkage awaits long reads.")
}

fwrite(data.table(
  quantity = c("alpha_sum_genomewide", "r_MOL", "r_PIP", "r_QUI", "pip_qui_gap",
               "pcor_depth_alphaQUI_given_bio05", "pcor_p",
               "pcor_depth_bio05_given_alphaQUI", "pcor_p_bio05"),
  value = c(s, r_MOL, r_PIP, r_QUI, r_QUI - r_PIP, pa$r, pa$p, pb$r, pb$p)),
  "results/chr1/Table_cyp6_depth_partial.csv")
say("")
say("wrote results/chr1/Table_cyp6_depth_partial.csv")
