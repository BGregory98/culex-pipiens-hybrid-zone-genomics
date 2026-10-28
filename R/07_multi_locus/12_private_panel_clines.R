# =============================================================================
# R/07_multi_locus/12_private_panel_clines.R
# =============================================================================
# Donor-private ancestry against the environment at the three candidate loci.
# A polarized axis (e.g. QUI vs MOL) leaves the third reference unconstrained,
# so a rise in the "QUI-type" allele cannot be apportioned between QUI arriving
# and MOL leaving. A private panel -- one reference >= 0.95 with BOTH others
# <= 0.05 -- can only be supplied by that lineage, so each lineage gets its own
# read-out. Within each interval, the per-pool mean frequency of each panel
# (alpha) is regressed on the locus predictor across the 24 field pools.
#
# Manuscript: Results "Cx. quinquefasciatus ancestry and sequencing read depth at
#             a CYP6 cluster rise with summer temperature" and "Cx.
#             quinquefasciatus ancestry is excluded at a developmental-quiescence
#             locus on the same gradient" (panel sizes, six coolest / six warmest
#             means, slopes); Supplementary Results "Donor-private panels give
#             opposing ancestry gradients at CpL1a and CpL2a"
#
# Inputs:  data/site_coords_with_spatial_data_allvars.csv   (R/00_setup/01_site_predictors.R)
#          deposit: data/allele_frequencies_v4.csv          (R/01_filtering/01_allele_frequency_table.R)
# Outputs: results/multilocus/Table_private_panel_alpha_per_pool.csv
#          results/chr1/Table_cpl1a_interval_all_snps.csv
#          results/chr2/Table_cpl2a_interval_all_snps.csv
#
# Runtime: under a minute (one pass over the allele-frequency table)
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")
suppressMessages(library(data.table))

AF_PATH <- deposit_file("data/allele_frequencies_v4.csv")
need_files("data/site_coords_with_spatial_data_allvars.csv")

HI <- 0.95; LO <- 0.05

# CpL2a is the 500 kb candidate core (13.5-14.0 Mb) in the private-panel analyses.
# Interval SNP tables are written only for the two narrow loci.
LOCI <- list(
  CpL1a = list(chr = CHR1, start = CPL1A_START, end = CPL1A_END,
               env = "bio05", cache = "results/chr1/Table_cpl1a_interval_all_snps.csv"),
  CpL2a = list(chr = CHR2, start = CPL2A_CORE_START,      end = CPL2A_CORE_END,
               env = "bio05", cache = "results/chr2/Table_cpl2a_interval_all_snps.csv"),
  # CpL3a spans 11.6 Mb / 245k SNPs; not written, it would be a ~35 MB table.
  CpL3a = list(chr = CHR3, start = CPL3A_START, end = CPL3A_END,
               env = "pis",   cache = NA)
)

hdr   <- names(fread(AF_PATH, nrows = 0, showProgress = FALSE))
pools <- setdiff(hdr, c("chromosome", "position", "cumulative_position"))
field <- setdiff(pools, REFS)

env <- fread("data/site_coords_with_spatial_data_allvars.csv", showProgress = FALSE)
setkey(env, site); env <- env[field]

say("reading %s (one pass, three intervals) ...", AF_PATH)
af <- fread(AF_PATH, select = c("chromosome", "position", pools), showProgress = FALSE)
say("matrix rows: %s", format(nrow(af), big.mark = ","))

out <- list(); alphas <- list()
for (nm in names(LOCI)) {
  L  <- LOCI[[nm]]
  iv <- af[chromosome == L$chr & position >= L$start & position <= L$end]
  x  <- as.numeric(env[[L$env]])

  q <- iv[[QUI]]; p <- iv[[PIP]]; m <- iv[[MOL]]
  ok <- is.finite(q) & is.finite(p) & is.finite(m)
  panels <- list(
    "QUI-private" = ok & q >= HI & p <= LO & m <= LO,
    "PIP-private" = ok & p >= HI & q <= LO & m <= LO,
    "MOL-private" = ok & m >= HI & q <= LO & p <= LO
  )
  F <- as.matrix(iv[, ..field])
  lo6 <- order(x)[1:6]; hi6 <- order(x, decreasing = TRUE)[1:6]

  say("")
  say("################ %s  (%.2f-%.2f Mb, predictor %s, %s SNPs) ################",
      nm, L$start/1e6, L$end/1e6, L$env, format(nrow(iv), big.mark = ","))
  say("%-12s %6s %9s %9s %11s %8s %11s", "panel", "n", "low env", "high env", "slope", "adjR2", "P")
  for (pn in names(panels)) {
    sel <- panels[[pn]]
    if (sum(sel) < 5) { say("%-12s %6d   too few SNPs to analyse", pn, sum(sel)); next }
    alpha <- colMeans(F[sel, , drop = FALSE], na.rm = TRUE)
    fit <- lm(alpha ~ x); s <- summary(fit)
    say("%-12s %6d %9.3f %9.3f %+11.5f %8.3f %11.3g",
        pn, sum(sel), mean(alpha[lo6]), mean(alpha[hi6]),
        coef(fit)[2], s$adj.r.squared, coef(s)[2, 4])
    alphas[[paste(nm, pn)]] <- data.table(locus = nm, panel = pn, pool = field,
                                          env = x, alpha = as.numeric(alpha))
  }
  # purity: mean reference frequencies over each panel's SNPs
  for (pn in names(panels)) {
    sel <- panels[[pn]]
    if (sum(sel) < 5) next
    say("   purity %-12s QUI %.3f  PIP %.3f  MOL %.3f", pn, mean(q[sel]), mean(p[sel]), mean(m[sel]))
  }
  if (!is.na(L$cache)) {
    dir.create(dirname(L$cache), showWarnings = FALSE, recursive = TRUE)
    fwrite(iv, L$cache)
    say("   wrote interval SNPs -> %s", L$cache)
  }
  rm(iv, F); invisible(gc())
}

res <- rbindlist(alphas)
fwrite(res, "results/multilocus/Table_private_panel_alpha_per_pool.csv")
say("")
say("wrote results/multilocus/Table_private_panel_alpha_per_pool.csv (%d rows)", nrow(res))
say("done.")
