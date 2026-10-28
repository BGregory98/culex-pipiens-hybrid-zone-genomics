# =============================================================================
# R/08_loci/09_cyp6_reference_pool_depth.R
# =============================================================================
# Cluster depth ratio at the CpL1a CYP6 paralogs in the three reference pools,
# against the field pools, from the callable-depth products re-aggregated with
# all 27 pools. Diagnoses reference-side inflation (collapsed paralogs in the
# assembly would lift all three references alike) against a lineage difference
# in copy number, fits field depth to the frequency of the CpL1a QUI-state
# haplotype and checks that line against the three references at their own
# haplotype frequencies, and rescales field ratios between the lowest and
# highest reference.
#
# Manuscript: Results "Cx. quinquefasciatus ancestry and sequencing read depth
#             at a CYP6 cluster rise with summer temperature" (0.999 / 1.514 /
#             1.579); Supplementary Results "Read depth in the reference pools,
#             and the haplotype that carries the extra copies"; Supplementary
#             Methods "Callable-position read-depth analysis"
#
# Inputs:  data/glmm_v4/cds_depth_per_gene_pool_all27.tsv          (cluster/06_depth/glmm/04b_aggregate_cds_depth_all27.slurm)
#          data/glmm_v4/genome_callable_depth_per_pool_all27.tsv   (cluster/06_depth/glmm/04b_aggregate_cds_depth_all27.slurm)
#          data/site_coords_with_spatial_data_allvars.csv          (tracked input; no producing script in R/)
#          data/sample_categories.csv                              (tracked input)
#          results/multilocus/Table_private_panel_alpha_per_pool.csv   (R/07_multi_locus/12_private_panel_clines.R)
#          data/focal_snp_results_cluster1.csv                     (R/07_multi_locus/06_allele_frequency_clines.R)
#          data/candidate_afs_cluster1.csv                         (R/07_multi_locus/06_allele_frequency_clines.R)
# Outputs: results/chr1/Table_cyp6_cluster_depth_with_references.csv
#          results/chr1/Table_cyp6_haplotype_depth.csv      (per pool: haplotype frequency,
#                                                            cluster ratio, fitted ratio)
#          results/chr1/Table_cyp6_haplotype_depth_fit.csv  (the field fit and the
#                                                            reference predictions)
#
# Runtime: < 1 min
# =============================================================================
#
# WHY THIS CONTROL
# ----------------
# The cluster sits above genome-average depth in every field pool (ratio 1.06 to
# 1.64) and the ratio rises with temperature. The artifact alternative is that
# reads not originating from these six loci pile onto them, which is live because
# the assembly is a Cx. pipiens f. pipiens individual, QUI is the most divergent
# lineage, and QUI ancestry rises with temperature -- so the pools with the
# highest apparent depth are also the pools whose reads match the reference
# least well. Field pools alone cannot separate those.
#
# HOW IT DISCRIMINATES
# --------------------
# PIP/ME1 is near-conspecific with the assembly and so carries almost no
# divergence penalty:
#
#   ratio_PIP near 1.0  -> the cluster is NOT intrinsically depth-inflated in
#                          this assembly. Elevation in the field pools is a
#                          property of the sample, and the excess must come from
#                          their non-PIP component. Supports the introgressed-
#                          duplication reading.
#   ratio_PIP well above 1 -> the cluster reads as depth-inflated even for a
#                          sample matching the assembly, i.e. collapsed paralogs
#                          baked into the reference. That would undercut the
#                          depth argument at its base, because "above
#                          genome-average" would be a fact about the denominator.
#
# QUI/FL1 alone is NOT decisive: QUI is both the most divergent from the
# assembly and the putative duplication donor, so a high QUI ratio is predicted
# by both hypotheses. What the three references jointly buy is a calibration --
# three near-pure endpoints at known divergence and known ancestry, against
# which the field pools can be checked.
#
# The callable set is defined by min-cov 10 across all 27 pools, so field and
# reference depths are measured at exactly the same positions. That is what makes
# the comparison clean, and it also means dropout is already excluded from both
# sides; this control speaks to copy number, not to dropout.
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")

suppressMessages(library(data.table))
# genome_callable_depth exceeds int32, so fread returns integer64. Without bit64
# loaded those values PRINT as garbage (2.4e-314). Only the printout is affected --
# total_depth and n_callable_pos are plain integers and the ratios below are
# computed from those -- but load bit64 when available so the log is readable.
if (requireNamespace("bit64", quietly = TRUE)) suppressMessages(library(bit64))

CDS_PATH   <- "data/glmm_v4/cds_depth_per_gene_pool_all27.tsv"
GW_PATH    <- "data/glmm_v4/genome_callable_depth_per_pool_all27.tsv"
ENV_FILE   <- "data/site_coords_with_spatial_data_allvars.csv"
ALPHA_FILE <- "results/multilocus/Table_private_panel_alpha_per_pool.csv"
FOCAL_FILE <- "data/focal_snp_results_cluster1.csv"
CAND_FILE  <- "data/candidate_afs_cluster1.csv"
OUT_PATH   <- "results/chr1/Table_cyp6_cluster_depth_with_references.csv"
HAP_OUT     <- "results/chr1/Table_cyp6_haplotype_depth.csv"
HAP_FIT_OUT <- "results/chr1/Table_cyp6_haplotype_depth_fit.csv"
need_files(CDS_PATH, GW_PATH, ENV_FILE, ALPHA_FILE, FOCAL_FILE, CAND_FILE)

REF <- c(QUI = QUI, PIP = PIP, MOL = MOL)

d   <- fread(CDS_PATH, showProgress = FALSE)
gw  <- fread(GW_PATH,  showProgress = FALSE)
env <- fread(ENV_FILE, showProgress = FALSE)
field <- field_names

pools <- sort(unique(d$pool))
say("pools in the depth table: %d (%d field, %d non-field)",
    length(pools), sum(pools %in% field), sum(!pools %in% field))
missing_ref <- setdiff(REF, pools)
if (length(missing_ref)) stop("reference pools absent from the table: ",
                              paste(missing_ref, collapse = ", "))

# ---- Normalise exactly as the per-paralog slopes (07_cyp6_copy_number.R) ----
# Depth per callable position, divided by the MEDIAN of that quantity across all
# genes in the same pool. The denominator is not free to choose: a depth-weighted
# mean gives materially different slopes and would not be comparable to the
# published values.
d[, per_pos := total_depth / n_callable_pos]
d[, ratio := per_pos / median(per_pos, na.rm = TRUE), by = pool]

present <- intersect(CYP6_GENES, unique(d$gene_id))
if (length(present) < length(CYP6_GENES))
  say("WARNING: only %d of the 6 focal paralogs present: %s",
      length(present), paste(present, collapse = ", "))

clu <- d[gene_id %in% present, .(cluster_ratio = mean(ratio)), by = pool]
setkey(clu, pool)

# =============================================================================
# 1. The three reference pools, and the PIP baseline
# =============================================================================
say("")
say("=== 1. cluster depth ratio in the three reference pools ===")
say("%-8s %-6s %14s %16s", "lineage", "pool", "cluster ratio", "gw callable depth")
for (l in names(REF)) {
  p <- REF[[l]]
  say("%-8s %-6s %14.3f %16s", l, p, clu[.(p), cluster_ratio],
      format(gw[pool == p, genome_callable_depth], big.mark = ","))
}
r_pip <- clu[.(REF[["PIP"]]), cluster_ratio]
r_qui <- clu[.(REF[["QUI"]]), cluster_ratio]
r_mol <- clu[.(REF[["MOL"]]), cluster_ratio]

fr <- clu[pool %in% field, cluster_ratio]
say("")
say("field pools for comparison: %.3f to %.3f (median %.3f)",
    min(fr), max(fr), median(fr))
say("")
say("--- diagnosis from the three-way pattern ---")
# Reference-side inflation from collapsed paralogs would raise ALL THREE
# references together, so it is the pattern across the three, not PIP's level by
# itself, that identifies the cause. All cases are enumerated.
ELEV <- 1.20   # a reference is "elevated" above this
SPREAD <- 0.25 # two references differ meaningfully by more than this
rr <- c(QUI = r_qui, PIP = r_pip, MOL = r_mol)
hi <- names(rr)[rr >= ELEV]; lo <- names(rr)[rr < ELEV]
say("QUI %.3f | PIP %.3f | MOL %.3f   (elevated: %s | at background: %s)",
    r_qui, r_pip, r_mol,
    if (length(hi)) paste(hi, collapse = ",") else "none",
    if (length(lo)) paste(lo, collapse = ",") else "none")
say("")
if (length(hi) == 0) {
  say("None of the references is elevated. The cluster is at background depth in")
  say("all three lineages, so the field-pool elevation is not lineage copy number.")
} else if (length(hi) == 3 && diff(range(rr)) < SPREAD) {
  say("All three references are elevated and close together (spread %.3f). That is",
      diff(range(rr)))
  say("the signature of reference-side inflation -- collapsed paralogs in the")
  say("assembly -- because a biological duplication would not be shared equally by")
  say("all three lineages. Absolute ratios should NOT be read as copy number.")
} else if (length(hi) == 1) {
  say("Only %s is elevated. That is the clean case: the high-copy state is specific", hi)
  say("to one lineage, and depth can attribute the field-pool elevation to it.")
} else {
  say("%s are elevated (spread %.3f) and %s sits at background. This is a genuine",
      paste(hi, collapse = " and "), diff(range(rr[hi])), paste(lo, collapse = ","))
  say("lineage difference in copy number rather than reference-side inflation --")
  say("inflation would lift all three. But it means depth CANNOT DISCRIMINATE the")
  say("elevated lineages from each other. A field pool gaining %s and a field pool",
      hi[1])
  say("gaining %s would show the same depth rise, so a depth cline supports 'losing",
      hi[2])
  say("%s' without identifying which of %s replaces it.",
      paste(lo, collapse = ","), paste(hi, collapse = "/"))
}
say("")
say("field median %.3f relative to each reference: QUI %+.3f | PIP %+.3f | MOL %+.3f",
    median(fr), median(fr) - r_qui, median(fr) - r_pip, median(fr) - r_mol)

# =============================================================================
# 2. Per-paralog detail, references beside the field range
# =============================================================================
say("")
say("=== 2. per-paralog ratios ===")
say("%-14s %8s %8s %8s %14s", "gene", "QUI", "PIP", "MOL", "field min-max")
for (g in present) {
  gi <- d[gene_id == g]
  setkey(gi, pool)
  fv <- gi[pool %in% field, ratio]
  say("%-14s %8.3f %8.3f %8.3f %6.2f-%-6.2f", g,
      gi[.(REF[["QUI"]]), ratio], gi[.(REF[["PIP"]]), ratio],
      gi[.(REF[["MOL"]]), ratio], min(fv), max(fv))
}

# =============================================================================
# 3. Field and reference pools against the CpL1a haplotype
# =============================================================================
# The 544 candidate SNPs nearly fixed between the QUI and MOL references
# (07_multi_locus/06_allele_frequency_clines.R) split CpL1a into two haplotype
# classes: the allele states QUI carries and those MOL carries. A pool's
# haplotype frequency is its mean frequency of the QUI-state allele over the 544
# SNPs, from read counts in the field pools and from the reference frequencies
# in the three references. If the extra cluster copies travel on the QUI-state
# haplotype, field cluster depth is linear in that frequency and the line fitted
# to the field pools passes through each reference at the reference's own
# frequency. Panel alpha cannot calibrate this: it is an index rather than a
# proportion, and PIP is not a single copy state at this locus (it carries the
# QUI-state haplotype at about 0.58).
say("")
say("=== 3. field and reference pools against the CpL1a haplotype ===")
fs  <- fread(FOCAL_FILE, showProgress = FALSE)
caf <- fread(CAND_FILE, showProgress = FALSE)[snp_id %in% fs$snp_id]
stopifnot(all(fs$higher_af %in% c(REF[["QUI"]], REF[["MOL"]])))
caf <- merge(caf, fs[, .(snp_id, q_alt = higher_af == REF[["QUI"]])], by = "snp_id")
caf[, q_state := fifelse(q_alt, frequency, 1 - frequency)]
hap <- caf[, .(hap_freq = mean(q_state), n_snps = .N), by = .(pool = sample)]
stopifnot(all(hap$n_snps == nrow(fs)), setequal(hap$pool, field))
ref_hap <- sapply(REF, function(p) {
  f <- fs[[paste0(p, "_freq")]]
  mean(fifelse(fs$higher_af == REF[["QUI"]], f, 1 - f))
})
say("%d SNPs; QUI-state haplotype frequency: field %.3f to %.3f; QUI %.3f, PIP %.3f, MOL %.3f",
    nrow(fs), min(hap$hap_freq), max(hap$hap_freq),
    ref_hap[["QUI"]], ref_hap[["PIP"]], ref_hap[["MOL"]])

al <- fread(ALPHA_FILE, showProgress = FALSE)
aq <- al[locus == "CpL1a" & panel == "QUI-private", .(pool, alpha_qui = alpha)]
fd <- merge(merge(clu[pool %in% field], hap[, .(pool, hap_freq)], by = "pool"),
            aq, by = "pool")
setkey(env, site)
fd[, bio05 := as.numeric(env[pool][["bio05"]])]
stopifnot(nrow(fd) == length(field))

f_h <- lm(cluster_ratio ~ hap_freq, data = fd)
sh  <- summary(f_h)
say("field: cluster_ratio = %.3f + %.3f x haplotype frequency (adj R2 %.3f, P %.3g)",
    coef(f_h)[1], coef(f_h)[2], sh$adj.r.squared, sh$coefficients[2, 4])
say("  at frequency 1 the line gives %.3f, %.2f times its value at frequency 0",
    sum(coef(f_h)), sum(coef(f_h)) / coef(f_h)[1])
pr <- predict(f_h, newdata = data.frame(hap_freq = ref_hap), interval = "prediction")
rownames(pr) <- names(REF)
for (l in names(REF))
  say("  %-3s at frequency %.3f: predicted %.3f (95%% PI %.3f-%.3f), observed %.3f",
      l, ref_hap[[l]], pr[l, "fit"], pr[l, "lwr"], pr[l, "upr"], rr[[l]])
f_hb <- lm(cluster_ratio ~ hap_freq + bio05, data = fd)
f_ha <- lm(cluster_ratio ~ hap_freq + alpha_qui, data = fd)
f_ab <- lm(hap_freq ~ bio05, data = fd)
say("with haplotype frequency in the model: bio05 P = %.3g, QUI-private alpha P = %.3g",
    summary(f_hb)$coefficients["bio05", 4], summary(f_ha)$coefficients["alpha_qui", 4])
say("r(haplotype frequency, QUI-private alpha) = %.3f; haplotype frequency ~ bio05 adj R2 %.3f",
    cor(fd$hap_freq, fd$alpha_qui), summary(f_ab)$adj.r.squared)

hap_out <- rbind(
  fd[, .(pool, pool_type = "field", hap_freq, cluster_ratio,
         predicted = unname(fitted(f_h)))],
  data.table(pool = unname(REF), pool_type = names(REF), hap_freq = unname(ref_hap),
             cluster_ratio = unname(rr[names(REF)]), predicted = unname(pr[names(REF), "fit"])))
fwrite(hap_out, HAP_OUT)
fwrite(data.table(
  quantity = c("n_snps", "intercept", "slope", "adj_r2", "p_slope", "ratio_at_fixation",
               "fold_at_fixation",
               paste0(rep(names(REF), each = 4), c("_hap_freq", "_predicted", "_pi_lo", "_pi_hi")),
               "p_bio05_given_hap", "p_alphaQUI_given_hap", "r_hap_alphaQUI"),
  value = c(nrow(fs), coef(f_h)[1], coef(f_h)[2], sh$adj.r.squared, sh$coefficients[2, 4],
            sum(coef(f_h)), sum(coef(f_h)) / coef(f_h)[1],
            as.vector(t(cbind(ref_hap[names(REF)], pr[names(REF), c("fit", "lwr", "upr")]))),
            summary(f_hb)$coefficients["bio05", 4],
            summary(f_ha)$coefficients["alpha_qui", 4], cor(fd$hap_freq, fd$alpha_qui))),
  HAP_FIT_OUT)
say("wrote %s and %s", HAP_OUT, HAP_FIT_OUT)

# implied duplication frequency under a two-endpoint PIP-to-QUI model, but only
# if that model is coherent (QUI must sit above PIP for the denominator to mean
# anything)
say("")
# The endpoints of this scale are the lowest and highest reference, not PIP and
# QUI by assumption. And the gap has to be wide enough to divide by: with PIP as
# the low endpoint, the QUI-PIP gap of 0.065 returns implied frequencies from
# -7.0 to +1.9. A scale whose denominator is a rounding error is not a scale.
lo_ref <- names(rr)[which.min(rr)]; hi_ref <- names(rr)[which.max(rr)]
r_lo <- min(rr); r_hi <- max(rr)
MIN_GAP <- 0.30
if (r_hi - r_lo > MIN_GAP) {
  fd[, implied_freq := (cluster_ratio - r_lo) / (r_hi - r_lo)]
  say("=== implied high-copy frequency, %s-to-%s scale (%.3f to %.3f) ===",
      lo_ref, hi_ref, r_lo, r_hi)
  say("Endpoints are the lowest and highest reference, not assumed to be PIP and QUI.")
  n_out <- sum(fd$implied_freq < -0.05 | fd$implied_freq > 1.05)
  if (n_out > 0)
    say("WARNING: %d of %d field pools fall outside [0,1]; the endpoints do not",
        n_out, nrow(fd))
  if (n_out > 0)
    say("bracket the field range, so read these as relative position, not frequency.")
  setorder(fd, bio05)
  for (i in seq_len(nrow(fd)))
    say("  %-12s bio05 %.3f  ratio %.3f  implied %.2f", fd$pool[i], fd$bio05[i],
        fd$cluster_ratio[i], fd$implied_freq[i])
  ff <- lm(implied_freq ~ bio05, data = fd)
  say("  implied frequency rises %+.3f per C (P = %.3g)",
      coef(ff)[2], summary(ff)$coefficients[2, 4])
  say("  NOTE: this rescales the SAME data and adds no independent evidence; it is")
  say("  reported because it is the quantity the copy-number claim is about.")
} else {
  say("=== two-endpoint scale not applicable ===")
  say("The reference spread is only %.3f (%s %.3f to %s %.3f), below the %.2f",
      r_hi - r_lo, lo_ref, r_lo, hi_ref, r_hi, MIN_GAP)
  say("minimum, so there is no denominator to build a frequency scale on.")
}
# The scale above is only meaningful between the two lineages it is drawn between.
# If the two elevated references are close to each other, the informative contrast
# is elevated-versus-background, and depth cannot say which elevated lineage a
# field pool is gaining.
if (length(hi) >= 2 && diff(range(rr[hi])) < SPREAD) {
  say("")
  say("--- what this scale can and cannot say ---")
  say("%s differ by only %.3f, so any scale drawn between them is noise. The",
      paste(hi, collapse = " and "), diff(range(rr[hi])))
  say("interpretable contrast is %s-versus-%s, and a field pool's position on it",
      paste(hi, collapse = "/"), paste(lo, collapse = ","))
  say("says how much %s-like sequence it has lost, NOT which lineage replaced it.",
      paste(lo, collapse = ","))
}

out <- merge(clu, gw[, .(pool, genome_callable_depth)], by = "pool", all.x = TRUE)
out[, pool_type := ifelse(pool %in% field, "field",
                   ifelse(pool %in% REF, names(REF)[match(pool, REF)], "other"))]
fwrite(out, OUT_PATH)
say("")
say("wrote %s", OUT_PATH)
say("")
say("Whatever this returns, it refines the mechanism behind one of three")
say("concordant signals at CpL1a. The locus-level finding rests on allele")
say("frequency and lineage-private ancestry, neither of which uses depth.")
