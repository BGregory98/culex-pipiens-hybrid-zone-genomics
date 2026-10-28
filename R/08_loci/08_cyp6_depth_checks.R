# =============================================================================
# R/08_loci/08_cyp6_depth_checks.R
# =============================================================================
# Tests whether the CpL1a depth cline could be a mapping artifact, using the
# 24 field pools:
#
#   1. Conservation. Cross-paralog mismapping within the cluster REDISTRIBUTES
#      reads among the six paralogs; it cannot raise the cluster total. If the
#      six rise together, within-cluster mismapping is excluded as the sole
#      cause; if they trade off, it is implicated.
#   2. Ancestry. If the depth cline is a by-product of QUI ancestry (the
#      mismapping route), alpha_QUI at the locus should absorb the bio05
#      signal. If depth carries information beyond ancestry, both survive a
#      joint fit.
#   3. Direction. Reference-mapping bias against QUI-like reads causes DROPOUT,
#      which lowers depth where QUI ancestry is high -- the opposite sign to
#      what is observed. Quantified as the sign of the depth-on-alpha_QUI slope.
#   4. Recruitment. Reads recruited into the cluster from elsewhere must be
#      lost by donor genes; the per-gene depth-temperature slopes of all genes
#      are checked for a heavier losing tail.
#
# Reference-pool depth at the cluster is in 09_cyp6_reference_pool_depth.R.
#
# Manuscript: Supplementary Results "Whether the CpL1a depth cline could be a
#             mismapping artifact"
#
# Inputs:  data/glmm_v4/cds_depth_per_gene_pool.tsv                (cluster/06_depth/glmm/04_aggregate_cds_depth.R)
#          data/site_coords_with_spatial_data_allvars.csv          (tracked input; no producing script in R/)
#          results/multilocus/Table_private_panel_alpha_per_pool.csv   (R/07_multi_locus/12_private_panel_clines.R)
# Outputs: results/chr1/Table_cyp6_cluster_depth_vs_ancestry.csv
#          results/chr1/Table_depth_slope_all_genes.csv
#
# Runtime: < 1 min
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")

suppressMessages(library(data.table))

DEPTH_FILE <- "data/glmm_v4/cds_depth_per_gene_pool.tsv"
ENV_FILE   <- "data/site_coords_with_spatial_data_allvars.csv"
ALPHA_FILE <- "results/multilocus/Table_private_panel_alpha_per_pool.csv"
need_files(DEPTH_FILE, ENV_FILE, ALPHA_FILE)

# ---- per-gene CDS depth, normalised by each pool's median gene ----
d   <- fread(DEPTH_FILE, showProgress = FALSE)
env <- fread(ENV_FILE, showProgress = FALSE)
al  <- fread(ALPHA_FILE, showProgress = FALSE)
aq <- al[locus == "CpL1a" & panel == "QUI-private", .(pool, alpha_qui = alpha)]

# Depth per callable position, divided by the MEDIAN of that quantity across all
# genes in the same pool. This is deliberately the same denominator as the
# per-paralog slopes of 07_cyp6_copy_number.R -- a depth-weighted mean denominator
# gives materially different slopes (e.g. CQUJHB009005 at +0.05 rather than +0.13),
# so using the mean here would contradict the reported values rather than test them.
d[, per_pos := total_depth / n_callable_pos]
d[, ratio := per_pos / median(per_pos, na.rm = TRUE), by = pool]

setkey(env, site)
pools <- sort(unique(d$pool))
x   <- as.numeric(env[pools][["bio05"]])
setkey(aq, pool)
aqv <- aq[pools][["alpha_qui"]]
stopifnot(length(x) == length(pools), !anyNA(x), !anyNA(aqv))

# =============================================================================
# 1. Conservation: do the six paralogs rise together, or trade off?
# =============================================================================
say("=== 1: do the six CYP6 paralogs rise together, or redistribute? ===")
say("Cross-paralog mismapping moves reads WITHIN the cluster. It changes the shares")
say("but cannot change the cluster total. Signs below are per degree C.")
say("")
say("%-16s %8s %8s %9s %9s %11s", "gene", "cool", "warm", "slope", "adjR2", "P")
M <- matrix(NA_real_, nrow = length(pools), ncol = length(CYP6_GENES),
            dimnames = list(pools, CYP6_GENES))
for (g in CYP6_GENES) {
  y <- d[gene_id == g][match(pools, pool)][["ratio"]]
  M[, g] <- y
  f <- lm(y ~ x); ss <- summary(f)
  say("%-16s %8.3f %8.3f %+9.4f %9.3f %11.3g", g,
      mean(y[order(x)[1:6]]), mean(y[order(x, decreasing = TRUE)[1:6]]),
      coef(f)[2], ss$adj.r.squared, ss$coefficients[2, 4])
}
nup <- sum(apply(M, 2, function(y) coef(lm(y ~ x))[2]) > 0)
say("")
say("paralogs with a POSITIVE depth-temperature slope: %d of %d", nup, length(CYP6_GENES))

tot <- rowSums(M)
f <- lm(tot ~ x); ss <- summary(f)
say("cluster TOTAL depth ratio: %.3f (coolest 6) -> %.3f (warmest 6), slope %+.4f/C, P = %.3g",
    mean(tot[order(x)[1:6]]), mean(tot[order(x, decreasing = TRUE)[1:6]]),
    coef(f)[2], ss$coefficients[2, 4])
say("")
say("pairwise correlations among the six paralogs' depth ratios across pools:")
cm <- cor(M)
say("  mean off-diagonal r = %+.3f ; range %+.3f to %+.3f ; negative pairs: %d of %d",
    mean(cm[upper.tri(cm)]), min(cm[upper.tri(cm)]), max(cm[upper.tri(cm)]),
    sum(cm[upper.tri(cm)] < 0), sum(upper.tri(cm)))
say("")
say("Reading: co-elevation with a rising TOTAL cannot be produced by redistribution")
say("among these six genes. Trading off (negative pairs, flat total) would be.")

# =============================================================================
# 2. Ancestry: does alpha_QUI absorb the temperature signal?
# =============================================================================
say("")
say("=== 2: depth against temperature, against QUI ancestry, and jointly ===")
say("%-30s %11s %11s %9s", "model", "P bio05", "P alpha_QUI", "adjR2")
m1 <- lm(tot ~ x); m2 <- lm(tot ~ aqv); m3 <- lm(tot ~ x + aqv)
say("%-30s %11.3g %11s %9.3f", "total ~ bio05",
    summary(m1)$coefficients[2, 4], "-", summary(m1)$adj.r.squared)
say("%-30s %11s %11.3g %9.3f", "total ~ alpha_QUI", "-",
    summary(m2)$coefficients[2, 4], summary(m2)$adj.r.squared)
say("%-30s %11.3g %11.3g %9.3f", "total ~ bio05 + alpha_QUI",
    summary(m3)$coefficients[2, 4], summary(m3)$coefficients[3, 4],
    summary(m3)$adj.r.squared)
say("")
say("r(bio05, alpha_QUI) = %+.3f -- the two predictors are collinear, so the joint",
    cor(x, aqv))
say("fit has limited power to separate them; report it as a direction check only.")

# =============================================================================
# 3. Direction: dropout vs inflation
# =============================================================================
say("")
say("=== 3: sign of the depth-on-QUI-ancestry relationship ===")
b <- coef(lm(tot ~ aqv))[2]
say("slope of cluster total depth on alpha_QUI: %+.3f", b)
say("Reference-mapping bias against QUI-like reads predicts DROPOUT (negative).")
say("Cross-paralog inflation predicts positive. Observed sign: %s",
    ifelse(b > 0, "POSITIVE -- not consistent with dropout", "negative"))
say("")
say("NOTE: a positive sign is consistent with BOTH real CNV on QUI-like haplotypes")
say("and with reads recruited into the cluster from OUTSIDE it. Reference-pool depth")
say("at the cluster (09_cyp6_reference_pool_depth.R) and section 4 address that.")

fwrite(data.table(pool = pools, bio05 = x, alpha_qui = aqv,
                  cluster_total_ratio = tot,
                  as.data.table(M)),
       "results/chr1/Table_cyp6_cluster_depth_vs_ancestry.csv")
say("")
say("wrote results/chr1/Table_cyp6_cluster_depth_vs_ancestry.csv")
# =============================================================================
# 4. Recruitment from outside the cluster: is there a donor set?
# =============================================================================
# The one route (1) cannot exclude is reads recruited into the cluster from
# homologous sequence elsewhere in the genome. No genome-wide P450 annotation
# is available, so a named homolog set cannot be tested. But the accounting constrains it without
# annotation: depth ratios are normalised per pool, so a gain of +0.93 units at the
# cluster must be matched by an equal loss elsewhere. Spread across 13,735 genes
# that is -7e-05 each and undetectable; concentrated in a handful of homologs it
# would produce a few genes with large NEGATIVE slopes, mirroring the focal six.
# So the question is whether the negative tail is heavier than the positive tail.
say("")
say("=== 4: is there a concentrated donor set losing depth? ===")
dw <- dcast(d, gene_id ~ pool, value.var = "ratio")
gid <- dw$gene_id
Mall <- as.matrix(dw[, ..pools])
keep <- rowSums(is.finite(Mall)) == length(pools)
gid <- gid[keep]; Mall <- Mall[keep, , drop = FALSE]
xc <- x - mean(x)
slopes <- as.vector((Mall %*% xc) / sum(xc^2))
isf <- gid %in% CYP6_GENES
say("genes with complete depth across all 24 pools: %s (focal six included)",
    format(length(gid), big.mark = ","))
say("")
say("focal six slopes: %s", paste(sprintf("%+.3f", slopes[isf]), collapse = " "))
say("  smallest focal slope: %+.4f per C", min(slopes[isf]))
say("")
say("in-genome null: where do the focal six rank among all genes?")
for (i in which(isf))
  say("  %-14s slope %+.4f   percentile %.2f%%", gid[i], slopes[i],
      100 * mean(slopes <= slopes[i]))
say("")
nf <- slopes[!isf]
thr <- min(slopes[isf])
say("tail symmetry, excluding the focal six (n = %s):", format(length(nf), big.mark = ","))
say("  genes with slope >= %+.3f (focal magnitude, gaining): %d", thr, sum(nf >= thr))
say("  genes with slope <= %+.3f (focal magnitude, losing):  %d", -thr, sum(nf <= -thr))
say("  summed depth-ratio units gained by that upper tail: %+.2f", sum(nf[nf >= thr]))
say("  summed depth-ratio units lost by that lower tail:   %+.2f", sum(nf[nf <= -thr]))
say("  mean %+.5f, median %+.5f, SD %.4f", mean(nf), median(nf), sd(nf))
say("  skewness of the non-focal slope distribution: %+.3f",
    mean((nf - mean(nf))^3) / sd(nf)^3)
say("")
say("Reading: a concentrated donor set would show a heavier lower tail than upper.")
say("A symmetric tail means the cluster's gain is not traceable to any small group")
say("of genes losing an equivalent amount, which is what recruitment would require.")
fwrite(data.table(gene_id = gid, focal = isf, slope_per_C = slopes),
       "results/chr1/Table_depth_slope_all_genes.csv")
say("")
say("wrote results/chr1/Table_depth_slope_all_genes.csv")
