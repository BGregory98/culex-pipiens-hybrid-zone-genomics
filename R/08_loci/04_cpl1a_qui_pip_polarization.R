# =============================================================================
# R/08_loci/04_cpl1a_qui_pip_polarization.R
# =============================================================================
# Re-polarizes the CpL1a candidate SNPs on QUI versus PIP instead of QUI versus
# MOL. Counts the SNPs informative on each contrast at a range of |dp|
# thresholds, reports PIP's frequency of the QUI-type allele at the QUI-MOL
# near-fixed SNPs, fits per-SNP OLS clines of polarized field-pool frequency on
# bio05 for each contrast, and summarizes the reference frequencies at the 13
# nonsynonymous CYP6 sites. The QUI-vs-MOL clines are refit with OLS on
# frequencies as a check on the binomial GLMs of 06_allele_frequency_clines.R
# (543 of 544 rising at P < 0.01).
#
# Manuscript: Results "Cx. quinquefasciatus ancestry and sequencing read depth
#             at a CYP6 cluster rise with summer temperature" (PIP carries the
#             QUI allele at 0.56-0.82); Supplementary Results "Reference
#             structure at the CpL1a candidate SNPs, and why the QUI-versus-PIP
#             contrast is unavailable"
#
# Inputs:  results/baypass/v4/baypass_outliers_clustered.csv   (R/03_baypass/01_baypass_outliers.R)
#          data/site_coords_with_spatial_data_allvars.csv      (tracked input; no producing script in R/)
#          results/chr1/Table_chr1_cds_consequences.csv        (03_cyp6_coding_consequences.R)
#          deposit: data/allele_frequencies_v4.csv
# Outputs: results/chr1/Table_cpl1a_candidate_pool_freqs.csv   (frequencies of the 654 candidates)
#
# Runtime: ~3 min (frequency table)
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")

suppressMessages(library(data.table))

CAND_FILE <- "results/baypass/v4/baypass_outliers_clustered.csv"
ENV_FILE  <- "data/site_coords_with_spatial_data_allvars.csv"
CONS_FILE <- "results/chr1/Table_chr1_cds_consequences.csv"
need_files(CAND_FILE, ENV_FILE, CONS_FILE)
AF_PATH   <- deposit_file("data/allele_frequencies_v4.csv")
FREQ_OUT  <- "results/chr1/Table_cpl1a_candidate_pool_freqs.csv"

# -----------------------------------------------------------------------------
# 1. Candidate SNPs in CpL1a
# -----------------------------------------------------------------------------
cand <- fread(CAND_FILE, select = c("chromosome", "position"), showProgress = FALSE)
cand <- unique(cand[chromosome == CHR1 & position > CPL1A_START & position < CPL1A_END])
say("CpL1a candidate SNPs: %d  (manuscript reports 654)", nrow(cand))

# -----------------------------------------------------------------------------
# 2. Allele frequencies for those SNPs
# -----------------------------------------------------------------------------
hdr   <- names(fread(AF_PATH, nrows = 0, showProgress = FALSE))
pools <- setdiff(hdr, c("chromosome", "position", "cumulative_position"))
field <- setdiff(pools, REFS)
say("pools in matrix: %d  (%d field + %d reference)", length(pools), length(field), length(REFS))
stopifnot(length(field) == 24, all(REFS %in% pools))

say("reading %s ...", AF_PATH)
af <- fread(AF_PATH, select = c("chromosome", "position", pools), showProgress = FALSE)
say("matrix rows: %s", format(nrow(af), big.mark = ","))
setkey(af, chromosome, position)
af <- af[cand, nomatch = 0]
# The candidate subset is written out so the checks below can be rerun without
# the full frequency table.
fwrite(af, FREQ_OUT)
say("wrote candidate frequencies to %s", FREQ_OUT)
say("candidate SNPs with allele frequencies: %d", nrow(af))

# -----------------------------------------------------------------------------
# 3. Environmental predictor
# -----------------------------------------------------------------------------
env <- fread(ENV_FILE, showProgress = FALSE)
env <- env[site %in% field]
say("field pools matched to bio05: %d of %d", nrow(env), length(field))
stopifnot(nrow(env) == length(field))
setkey(env, site)
env   <- env[field]                 # order to match the frequency columns
bio05 <- as.numeric(env$bio05)
say("bio05 range: %.2f to %.2f C (SD %.3f)", min(bio05), max(bio05), sd(bio05))

# -----------------------------------------------------------------------------
# 4. Reference frequencies and the two contrasts
# -----------------------------------------------------------------------------
fq  <- as.matrix(af[, ..field])     # SNPs x 24 field pools
qui <- as.numeric(af[[QUI]]); pip <- as.numeric(af[[PIP]]); mol <- as.numeric(af[[MOL]])

ok  <- is.finite(qui) & is.finite(pip) & is.finite(mol) &
       rowSums(!is.finite(fq)) == 0
say("SNPs with complete reference and field data: %d of %d", sum(ok), nrow(af))
fq <- fq[ok, , drop = FALSE]; qui <- qui[ok]; pip <- pip[ok]; mol <- mol[ok]

d_QM <- abs(qui - mol)
d_QP <- abs(qui - pip)

say("")
say("=== informative SNPs by contrast and threshold ===")
say("%-10s %10s %10s", "threshold", "QUI-MOL", "QUI-PIP")
for (th in c(0.95, 0.90, 0.80, 0.70, 0.50, 0.30)) {
  say("%-10.2f %10d %10d", th, sum(d_QM >= th), sum(d_QP >= th))
}

# -----------------------------------------------------------------------------
# 5. PIP's frequency of the "QUI-diagnostic" allele
# -----------------------------------------------------------------------------
sel_QM <- d_QM >= 0.95
# orient so the reported allele is the one more common in QUI than in MOL
sgn    <- ifelse(qui >= mol, 1, -1)
pip_as_qui_type <- ifelse(sgn > 0, pip, 1 - pip)
mol_as_qui_type <- ifelse(sgn > 0, mol, 1 - mol)
qui_as_qui_type <- ifelse(sgn > 0, qui, 1 - qui)

say("")
say("=== at the %d QUI-vs-MOL near-fixed SNPs, frequency of the QUI-type allele ===", sum(sel_QM))
q <- function(v) quantile(v, c(0, .25, .5, .75, 1), na.rm = TRUE)
say("  in QUI : min %.3f  q1 %.3f  median %.3f  q3 %.3f  max %.3f", q(qui_as_qui_type[sel_QM])[1], q(qui_as_qui_type[sel_QM])[2], q(qui_as_qui_type[sel_QM])[3], q(qui_as_qui_type[sel_QM])[4], q(qui_as_qui_type[sel_QM])[5])
say("  in PIP : min %.3f  q1 %.3f  median %.3f  q3 %.3f  max %.3f", q(pip_as_qui_type[sel_QM])[1], q(pip_as_qui_type[sel_QM])[2], q(pip_as_qui_type[sel_QM])[3], q(pip_as_qui_type[sel_QM])[4], q(pip_as_qui_type[sel_QM])[5])
say("  in MOL : min %.3f  q1 %.3f  median %.3f  q3 %.3f  max %.3f", q(mol_as_qui_type[sel_QM])[1], q(mol_as_qui_type[sel_QM])[2], q(mol_as_qui_type[sel_QM])[3], q(mol_as_qui_type[sel_QM])[4], q(mol_as_qui_type[sel_QM])[5])
say("  -> PIP carries the QUI-type allele at median %.3f; if this is high, these", median(pip_as_qui_type[sel_QM]))
say("     are MOL-absent alleles shared by QUI and PIP, not QUI-vs-PIP markers.")

# -----------------------------------------------------------------------------
# 6. Vectorised per-SNP OLS of polarised frequency on bio05
# -----------------------------------------------------------------------------
cline <- function(keep, orient_sgn, label) {
  Y  <- fq[keep, , drop = FALSE]
  s  <- orient_sgn[keep]
  # Polarise so high = target-type allele. Row-wise flip by subsetting, NOT
  # ifelse(): ifelse() returns an object shaped like its *test*, so passing a
  # length-nrow vector against a matrix silently collapses Y to a vector.
  flip <- which(s < 0)
  if (length(flip)) Y[flip, ] <- 1 - Y[flip, , drop = FALSE]
  n  <- ncol(Y)
  xc <- bio05 - mean(bio05)
  Sxx <- sum(xc^2)
  ybar  <- rowMeans(Y)
  slope <- as.numeric((Y - ybar) %*% xc) / Sxx
  fitted <- outer(slope, xc) + ybar
  rss    <- rowSums((Y - fitted)^2)
  se     <- sqrt((rss / (n - 2)) / Sxx)
  tstat  <- slope / se
  pval   <- 2 * pt(-abs(tstat), df = n - 2)
  pval[!is.finite(pval)] <- 1

  say("")
  say("=== %s: %d SNPs ===", label, nrow(Y))
  say("  significant at P<0.05 : %d rising, %d falling", sum(pval < .05 & slope > 0), sum(pval < .05 & slope < 0))
  say("  significant at P<0.01 : %d rising, %d falling", sum(pval < .01 & slope > 0), sum(pval < .01 & slope < 0))
  say("  median slope          : %+.4f per C   (all SNPs)", median(slope))
  say("  median predicted freq : %.3f at coolest site, %.3f at warmest",
      median(ybar + slope * (min(bio05) - mean(bio05))),
      median(ybar + slope * (max(bio05) - mean(bio05))))
  invisible(list(slope = slope, pval = pval, Y = Y))
}

r_QM <- cline(sel_QM, sgn, "QUI-vs-MOL polarisation, |dp| >= 0.95 (published contrast)")
say("  [binomial GLMs on read counts (06_allele_frequency_clines.R): 543 of 544 rising at P < 0.01, none falling]")

# Where do the three references sit on the QUI-vs-MOL polarised scale?
say("")
say("=== reference positions on the QUI-vs-MOL polarised scale (%d SNPs) ===", sum(sel_QM))
say("  QUI median %.3f | PIP median %.3f | MOL median %.3f",
    median(qui_as_qui_type[sel_QM]), median(pip_as_qui_type[sel_QM]), median(mol_as_qui_type[sel_QM]))

for (th in c(0.95, 0.70)) {
  sel <- d_QP >= th
  if (sum(sel) >= 3) {
    cline(sel, ifelse(qui >= pip, 1, -1),
          sprintf("QUI-vs-PIP polarisation, |dp| >= %.2f", th))
  } else {
    say("")
    say("=== QUI-vs-PIP, |dp| >= %.2f: only %d informative SNPs -- not analysable ===", th, sum(sel))
  }
}

# -----------------------------------------------------------------------------
# 7. The 13 nonsynonymous sites (Table S6). This table already carries the
#    three reference frequencies.
# -----------------------------------------------------------------------------
cons <- fread(CONS_FILE, showProgress = FALSE)
ns   <- cons[consequence == "nonsynonymous"]
say("")
say("=== Table S6: %d nonsynonymous candidate SNPs in CpL1a ===", nrow(ns))
# orient to the allele more common in QUI than in MOL, matching the paper's axis
s_ns  <- ifelse(ns$FL1_altfreq >= ns$CAL1_altfreq, 1, -1)
pip_ns <- ifelse(s_ns > 0, ns$ME1_altfreq, 1 - ns$ME1_altfreq)
qui_ns <- ifelse(s_ns > 0, ns$FL1_altfreq, 1 - ns$FL1_altfreq)
mol_ns <- ifelse(s_ns > 0, ns$CAL1_altfreq, 1 - ns$CAL1_altfreq)
say("  QUI-type allele in PIP: min %.2f  median %.2f  max %.2f",
    min(pip_ns), median(pip_ns), max(pip_ns))
say("  QUI median %.2f | MOL median %.2f", median(qui_ns), median(mol_ns))
say("  SNPs where PIP is nearer QUI than MOL: %d of %d",
    sum(abs(pip_ns - qui_ns) < abs(pip_ns - mol_ns)), nrow(ns))
say("  QUI-vs-PIP |dp| >= 0.95 among these: %d ; >= 0.70: %d",
    sum(abs(ns$FL1_altfreq - ns$ME1_altfreq) >= 0.95),
    sum(abs(ns$FL1_altfreq - ns$ME1_altfreq) >= 0.70))
say("  genes carrying them: %s", paste(sort(unique(ns$gene_id)), collapse = ", "))

say("")
say("done.")
