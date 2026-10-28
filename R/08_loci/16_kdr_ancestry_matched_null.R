# =============================================================================
# R/08_loci/16_kdr_ancestry_matched_null.R
# =============================================================================
# Tests the correlation of kdr L1014F frequency with each pool's genome-wide
# donor ancestry (and with CpL1a CYP6 cluster depth) against SNPs matched to
# kdr on differentiation. kdr is one site, so the null is SNP-level: kdr is a
# major XtXst outlier (128.16), and any SNP with strong among-population
# structure correlates with any ancestry axis more than an average SNP does.
# Comparison SNPs are matched on
#   XtXst   -- differentiation, the ascertainment concern
#   M_P     -- BayPass posterior mean allele frequency, because a SNP near
#              fixation has less room to correlate with anything
# at three tolerances, and the same Spearman correlation is computed for focal
# and null alike. A random genome-wide SNP sample then measures how ancestry
# tracking rises with XtXst.
#
# Manuscript: Supplementary Results "A characterized resistance allele outside
#             the candidate set" (null median |rho| 0.574, P = 0.72; MOL and
#             QUI P = 0.88 and 0.97; random-SNP |rho| 0.165, top XtXst decile
#             0.323; kdr against CYP6 depth raw P = 0.89, partial P = 0.053);
#             Supplementary Methods "Targeted analysis of insecticide-resistance
#             target sites"
#
# Inputs:  results/multilocus/per_pool_ancestry_per_pool.csv    (R/07_multi_locus/01_per_pool_ancestry.R)
#          results/resistance/Table_kdr_L1014F_by_site.csv      (14_kdr_frequency.R)
#          data/glmm_v4/cds_depth_per_gene_pool_all27.tsv       (cluster/06_depth/glmm/04b_aggregate_cds_depth_all27.slurm)
#          deposit: results/baypass/v4/baypass_core_envsubset_summary_pi_xtx.out
#          deposit: data/allele_frequencies_v4.csv
# Outputs: results/resistance/Table_kdr_ancestry_matched_null.csv
#          results/resistance/Table_xtxst_vs_ancestry_tracking.csv
#
# Runtime: ~5 min (BayPass summary and frequency table)
# =============================================================================
#
# WHY MATCHING ON XtXst IS NECESSARY AND NOT CIRCULAR. XtXst is already
# Omega-corrected -- it is a Mahalanobis-type distance under the population
# covariance -- so the reasonable objection is that it has conditioned on
# population structure and cannot be smuggling ancestry back in. Section (3)
# measures this directly: a random SNP tracks ancestry at about the chance
# level for n = 24, and ancestry tracking rises monotonically with XtXst. The
# reason is that XtXst is a MAGNITUDE, not a direction: it says a SNP deviates
# from the expected covariance more than most, not that it deviates along an
# axis orthogonal to ancestry. Omega is moreover estimated across all 27 pools,
# where the three reference colonies dominate the covariance, so it is
# calibrated on reference-versus-field separation and barely penalises
# deviation along the much finer ancestry gradient among the 24 field pools.
# Conditioning on XtXst therefore selects SNPs whose deviation is large, and
# among the field pools the axis with room to deviate is largely ancestry.
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")

suppressPackageStartupMessages({
  library(data.table)
})

ANC_FILE   <- 'results/multilocus/per_pool_ancestry_per_pool.csv'
KDR_FILE   <- 'results/resistance/Table_kdr_L1014F_by_site.csv'
DEPTH_FILE <- 'data/glmm_v4/cds_depth_per_gene_pool_all27.tsv'
need_files(ANC_FILE, KDR_FILE, DEPTH_FILE)
XTX_PATH <- deposit_file('results/baypass/v4/baypass_core_envsubset_summary_pi_xtx.out')
AF_PATH  <- deposit_file('data/allele_frequencies_v4.csv')

# BayPass marker of kdr L1014F. MRK is the SNP's row in allele_frequencies_v4.csv:
# the two tables have the same SNPs in the same order.
KDR_MRK  <- 7016575L
KDR_POS  <- 191086281L
RAND_FRAC <- 0.0011     # fraction of all SNPs drawn for the random genome-wide sample
set.seed(42)

# ---- per-pool ancestry axes -------------------------------------------------

anc <- fread(ANC_FILE)
bg  <- dcast(anc[panel %in% c('ME1_private', 'CAL1_private', 'FL1_private')],
             pool ~ panel, value.var = 'alpha_genome')
setnames(bg, c('ME1_private', 'CAL1_private', 'FL1_private'),
             c('alpha_PIP', 'alpha_MOL', 'alpha_QUI'), skip_absent = TRUE)
bg <- bg[pool %in% field_names]
setorder(bg, pool)
POOLS <- bg$pool
stopifnot(length(POOLS) == 24L)

# ---- matched and random SNP sets ---------------------------------------------

xtx <- fread(XTX_PATH, select = c('MRK', 'M_P', 'XtXst'), showProgress = FALSE)
stopifnot(identical(xtx$MRK, seq_len(nrow(xtx))))
say('kdr: XtXst %.2f, M_P %.3f', xtx$XtXst[KDR_MRK], xtx$M_P[KDR_MRK])

# Tolerances around kdr's XtXst (128.16) and posterior mean frequency (0.606),
# rounded: XtXst +-10% = 115.35-140.98, +-25% = 96.1-160.2; M_P +-0.05 and +-0.10.
sets <- list(
  `XtXst +-10%, M_P +-0.05` =
    xtx[XtXst >= 115.35 & XtXst <= 140.98 & M_P >= 0.556 & M_P <= 0.656, MRK],
  `XtXst +-25%, M_P +-0.10` =
    xtx[XtXst >= 96.1 & XtXst <= 160.2 & M_P >= 0.506 & M_P <= 0.706, MRK],
  `XtXst +-10%, freq unmatched` =
    xtx[XtXst >= 115.35 & XtXst <= 140.98, MRK])
for (sn in names(sets)) say('  %-28s %d SNPs', sn, length(sets[[sn]]))

# Random genome-wide sample: each SNP kept with probability RAND_FRAC.
ra <- xtx[which(runif(nrow(xtx)) < RAND_FRAC)]
say('random genome-wide sample: %d SNPs', nrow(ra))

# Per-pool frequencies for the matched and random SNPs, by row.
af <- fread(AF_PATH, select = c('chromosome', 'position', POOLS), showProgress = FALSE)
stopifnot(nrow(af) == nrow(xtx),
          af$chromosome[KDR_MRK] == CHR2, af$position[KDR_MRK] == KDR_POS)
all_mrk <- sort(unique(unlist(sets)))
stopifnot(KDR_MRK %in% all_mrk)
F24 <- as.matrix(af[all_mrk, ..POOLS])
rownames(F24) <- all_mrk
RF  <- as.matrix(af[ra$MRK, ..POOLS])
rm(af, xtx); invisible(gc())

# Spearman of a SNP's 24 field-pool frequencies against an ancestry axis.
rho_vs <- function(freqmat, y) {
  ry <- rank(y)
  apply(freqmat, 1, function(x) {
    if (!all(is.finite(x)) || sd(x) == 0) return(NA_real_)
    cor(rank(x), ry)
  })
}

# The CYP6 cluster depth vector, so the SNP-side null can be run on it too.
# 15_kdr_cpl1a_matched_null.R asks "given kdr, is the CYP6 cluster unusual among
# gene clusters?". The symmetric question -- "given CYP6, is kdr unusual among
# equally differentiated SNPs?" -- is asked here, because if any high-XtXst SNP
# tracks CYP6 depth equally well then the pairing is structure, not resistance
# biology.
dp <- fread(DEPTH_FILE)[n_callable_pos > 0]
dp[, per_pos := total_depth / n_callable_pos]
dp[, ratio   := per_pos / median(per_pos, na.rm = TRUE), by = pool]
cyp <- dp[gene_id %in% CYP6_GENES, .(cyp6 = mean(ratio)), by = pool]
kd  <- fread(KDR_FILE)[, .(pool = site, bio05)]
cv  <- Reduce(function(x, y) merge(x, y, by = 'pool'), list(bg, cyp, kd))
setorder(cv, pool)
stopifnot(identical(cv$pool, POOLS))
Z   <- as.matrix(cv[, .(bio05, alpha_PIP)])
rc  <- resid(lm(cv$cyp6 ~ Z))

partial_vs <- function(freqmat) apply(freqmat, 1, function(x) {
  if (!all(is.finite(x)) || sd(x) == 0) return(NA_real_)
  cor(resid(lm(x ~ Z)), rc)
})

axes <- list(alpha_PIP = rho_vs(F24, bg$alpha_PIP),
             alpha_MOL = rho_vs(F24, bg$alpha_MOL),
             alpha_QUI = rho_vs(F24, bg$alpha_QUI),
             `CYP6 depth (raw)`     = rho_vs(F24, cv$cyp6),
             `CYP6 depth (partial)` = partial_vs(F24))

res <- rbindlist(lapply(names(axes), function(ax) {
  r_all <- axes[[ax]]
  focal <- r_all[as.character(KDR_MRK)]
  rbindlist(lapply(names(sets), function(sn) {
    mrk  <- sets[[sn]]
    nullr <- r_all[as.character(setdiff(mrk, KDR_MRK))]
    nullr <- nullr[is.finite(nullr)]
    data.table(axis = ax, null_set = sn, n_null = length(nullr),
               focal_rho = round(focal, 3),
               null_median_abs = round(median(abs(nullr)), 3),
               emp_P = round((sum(abs(nullr) >= abs(focal)) + 1) / (length(nullr) + 1), 4))
  }))
}))

cat('\n=== kdr L1014F vs genome-wide ancestry, against XtXst/frequency-matched SNPs ===\n')
print(res, row.names = FALSE)
fwrite(res, 'results/resistance/Table_kdr_ancestry_matched_null.csv')

# Recomputing the focal from the frequency file should reproduce the value
# 14_kdr_frequency.R gets from the raw allele counts.
cat(sprintf('\nfocal rho vs alpha_PIP from frequency file: %+.3f  (14_kdr_frequency.R, from counts: -0.514)\n',
            rho_vs(F24, bg$alpha_PIP)[as.character(KDR_MRK)]))

# ---- (3) is the XtXst matching doing real work? ------------------------------
# The justification for the null design, measured rather than argued. See the
# header note.

ra[, absrho := abs(rho_vs(RF, bg$alpha_PIP))]
ra <- ra[is.finite(absrho)]
ra[, dec := cut(XtXst, quantile(XtXst, 0:10 / 10), include.lowest = TRUE, labels = FALSE)]
dtab <- ra[, .(n = .N, XtXst_median = round(median(XtXst), 1),
               median_abs_rho = round(median(absrho), 3)), by = dec][order(dec)]
cat(sprintf('\n=== does XtXst predict ancestry tracking? (%d random SNPs) ===\n', nrow(ra)))
cat(sprintf('median |rho| vs alpha_PIP, random SNP: %.3f (independence expectation %.3f)\n',
            median(ra$absrho), 0.8 / sqrt(length(POOLS))))
print(dtab, row.names = FALSE)
cat(sprintf('Spearman(XtXst, |rho|) = %+.3f ; kdr XtXst 128.16 = %.2f percentile\n',
            cor(ra$XtXst, ra$absrho, method = 'spearman'),
            100 * mean(ra$XtXst < 128.16)))
fwrite(dtab, 'results/resistance/Table_xtxst_vs_ancestry_tracking.csv')

cat('\nWrote results/resistance/Table_kdr_ancestry_matched_null.csv\n')
cat('Wrote results/resistance/Table_xtxst_vs_ancestry_tracking.csv\n')
