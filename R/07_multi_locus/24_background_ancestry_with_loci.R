# =============================================================================
# R/07_multi_locus/24_background_ancestry_with_loci.R
# =============================================================================
# Whether excluding the three candidate intervals shapes the genome-wide donor
# ancestry trends. Genome-wide alpha (alpha_genome) is measured at
# donor-private panel SNPs outside CpL1a, CpL2a and CpL3a; here per-pool alpha
# is computed twice, at panel SNPs outside the intervals (which must reproduce
# alpha_genome) and at all panel SNPs, and the models of
# 22_background_ancestry_vs_environment.R (alpha on bio05, on impervious
# surface, on both, and with a metropolitan-area term) are refitted on each.
#
# Panels use the definition of 12_private_panel_clines.R: one reference >= 0.95
# with both others <= 0.05. alpha_genome comes from the genome-wide null, which
# applies it to exact frequencies from read counts, AD / (RD + AD). Rounding to
# two decimals moves SNPs across the cutoffs (0.946 becomes 0.95), so the
# allele-frequency table alone does not reproduce alpha_genome. The table
# therefore selects the candidate rows (rounding can add SNPs to a panel but
# never remove one), and the count table supplies the read counts at those rows.
#
# Manuscript: Supplementary Results "Genome-wide molestus-like ancestry follows
#             impervious surface independently of position" (genome-wide
#             ancestry with and without the candidate intervals)
#
# Inputs:  deposit: data/allele_frequencies_v4.csv           (R/01_filtering/01_allele_frequency_table.R)
#          deposit: data/allele_counts_v4.tsv                (cluster/04_baypass/vcf_to_allele_counts.sh)
#          results/multilocus/genome_wide_panel_vs_env_per_pool.csv     (22_background_ancestry_vs_environment.R)
#          results/multilocus/Table_genome_wide_panel_vs_env.csv        (22_background_ancestry_vs_environment.R)
#          results/multilocus/Table_genome_wide_panel_vs_pis_robustness.csv (22_background_ancestry_vs_environment.R)
# Outputs: results/multilocus/Table_genome_wide_panel_loci_sensitivity.csv
#            (one row per donor x model: slope, SE, P, adjusted R2 outside the
#             intervals and with them, panel sizes, r between the two alphas)
#          results/multilocus/genome_wide_panel_loci_sensitivity_per_pool.csv
#            (one row per donor x pool: SNPs used and alpha with and without
#             the intervals)
#
# Runtime: about a minute (one pass over each deposit table; the count table is
#          read in chunks of two million rows)
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")
suppressMessages(library(data.table))

AF_PATH <- deposit_file("data/allele_frequencies_v4.csv")
CT_PATH <- deposit_file("data/allele_counts_v4.tsv")
PER_POOL_22 <- "results/multilocus/genome_wide_panel_vs_env_per_pool.csv"
MODELS_22   <- "results/multilocus/Table_genome_wide_panel_vs_env.csv"
ROBUST_22   <- "results/multilocus/Table_genome_wide_panel_vs_pis_robustness.csv"
need_files(PER_POOL_22, MODELS_22, ROBUST_22)

OUT_TABLE    <- "results/multilocus/Table_genome_wide_panel_loci_sensitivity.csv"
OUT_PER_POOL <- "results/multilocus/genome_wide_panel_loci_sensitivity_per_pool.csv"

HI <- 0.95; LO <- 0.05   # panel cutoffs (12_private_panel_clines.R)
MIN_DEPTH <- 10          # as in the genome-wide null (every pool reaches it at every panel SNP)
CHUNK <- 2e6

# ---- candidate panel rows from the allele-frequency table -------------------
hdr   <- names(fread(AF_PATH, nrows = 0, showProgress = FALSE))
pools <- setdiff(hdr, c("chromosome", "position", "cumulative_position"))
field <- setdiff(pools, REFS)
stopifnot(length(field) == 24, setequal(field, field_names))

say("reading %s ...", AF_PATH)
af <- fread(AF_PATH, select = c("chromosome", "position", REFS), showProgress = FALSE)
n_rows <- nrow(af)
q <- af[[QUI]]; p <- af[[PIP]]; m <- af[[MOL]]
rows <- which((q >= HI & p <= LO & m <= LO) |
              (p >= HI & q <= LO & m <= LO) |
              (m >= HI & q <= LO & p <= LO))
key <- af[rows, .(chromosome, position)]
rm(af, q, p, m); invisible(gc())
say("  %s SNPs; %s candidate panel rows at two-decimal frequencies",
    format(n_rows, big.mark = ","), format(length(rows), big.mark = ","))

# ---- read counts at those rows (same SNPs, same order, no header) ------------
say("reading %s in chunks ...", CT_PATH)
got <- list()
for (s in seq(0, n_rows - 1, by = CHUNK)) {
  ch <- fread(CT_PATH, sep = "\t", header = FALSE, skip = s, nrows = min(CHUNK, n_rows - s),
              colClasses = "character", showProgress = FALSE)
  got[[length(got) + 1]] <- ch[rows[rows > s & rows <= s + nrow(ch)] - s]
  rm(ch); invisible(gc())
}
ct <- rbindlist(got); rm(got)
stopifnot(nrow(ct) == length(rows), ncol(ct) == 4 + length(pools),
          identical(ct$V1, key$chromosome), identical(as.integer(ct$V2), key$position))
rd_ad <- as.matrix(ct[, -(1:4)]); colnames(rd_ad) <- pools; rm(ct)
rd  <- matrix(as.numeric(sub(":.*", "", rd_ad)), nrow(rd_ad), dimnames = list(NULL, pools))
alt <- matrix(as.numeric(sub(".*:", "", rd_ad)), nrow(rd_ad), dimnames = list(NULL, pools))
dp  <- rd + alt
freq <- alt / dp
freq[dp < MIN_DEPTH] <- NA_real_

# ---- panels and the candidate-interval mask ---------------------------------
fq <- freq[, QUI]; fp <- freq[, PIP]; fm <- freq[, MOL]
panels <- list(QUI = fq >= HI & fp <= LO & fm <= LO,
               PIP = fp >= HI & fq <= LO & fm <= LO,
               MOL = fm >= HI & fq <= LO & fp <= LO)
panels <- lapply(panels, function(v) { v[is.na(v)] <- FALSE; v })
in_loci <- Reduce(`|`, lapply(LOCI, function(L)
  key$chromosome == L$chr & key$position >= L$start & key$position <= L$end))

# ---- per-pool alpha, with and without the intervals --------------------------
pp22 <- fread(PER_POOL_22)
per_pool <- rbindlist(lapply(names(panels), function(dn) {
  out_sel <- panels[[dn]] & !in_loci
  all_sel <- panels[[dn]]
  F_out <- freq[out_sel, field, drop = FALSE]
  F_all <- freq[all_sel, field, drop = FALSE]
  data.table(donor = dn, pool = field,
             n_snps_outside   = colSums(!is.na(F_out)),
             n_snps_with_loci = colSums(!is.na(F_all)),
             alpha_outside    = colMeans(F_out, na.rm = TRUE),
             alpha_with_loci  = colMeans(F_all, na.rm = TRUE))
}))
per_pool[, alpha_diff := alpha_with_loci - alpha_outside]
# predictors and metropolitan area exactly as 22 used them
per_pool <- merge(per_pool, pp22[, .(donor, pool, panel, alpha_genome, bio05, pis, metro)],
                  by = c("donor", "pool"))
stopifnot(nrow(per_pool) == 72)

# The outside version must be alpha_genome.
dev <- max(abs(per_pool$alpha_outside - per_pool$alpha_genome))
say("outside-interval alpha vs alpha_genome (22): max |difference| = %.2g", dev)
if (dev > 1e-10) stop("outside-interval alpha does not reproduce alpha_genome", call. = FALSE)
per_pool[, alpha_genome := NULL]
setcolorder(per_pool, c("donor", "panel", "pool"))
setorder(per_pool, donor, bio05)
fwrite(per_pool, OUT_PER_POOL)

# ---- refit the models on both versions -------------------------------------
fit_row <- function(dd, y, formula) {
  dd <- copy(dd); dd[, alpha := get(y)]
  s  <- summary(lm(formula, data = dd))
  co <- coef(s)
  data.table(focal_term = rownames(co)[2], slope = co[2, 1], se = co[2, 2],
             p = co[2, 4], adj_r2 = s$adj.r.squared,
             other_terms = paste(rownames(co)[-(1:2)], collapse = ";"),
             other_p = paste(signif(co[-(1:2), 4], 3), collapse = ";"))
}
MODELS <- list(
  "alpha ~ bio05"               = alpha ~ bio05,
  "alpha ~ pis"                 = alpha ~ pis,
  "alpha ~ bio05 + pis"         = alpha ~ bio05 + pis,
  "alpha ~ pis + bio05"         = alpha ~ pis + bio05,
  "alpha ~ bio05 + metro"       = alpha ~ bio05 + metro,
  "alpha ~ pis + bio05 + metro" = alpha ~ pis + bio05 + metro)

tab <- rbindlist(lapply(c("QUI", "PIP", "MOL"), function(dn) {
  dd <- per_pool[donor == dn]
  rbindlist(lapply(names(MODELS), function(k) {
    a <- fit_row(dd, "alpha_outside",   MODELS[[k]])
    b <- fit_row(dd, "alpha_with_loci", MODELS[[k]])
    data.table(donor = dn, model = k, focal_term = a$focal_term, n_pools = nrow(dd),
               slope_outside = a$slope, se_outside = a$se, p_outside = a$p,
               adj_r2_outside = a$adj_r2,
               slope_with_loci = b$slope, se_with_loci = b$se, p_with_loci = b$p,
               adj_r2_with_loci = b$adj_r2,
               other_terms = a$other_terms, other_p_outside = a$other_p,
               other_p_with_loci = b$other_p,
               n_panel_outside = sum(panels[[dn]] & !in_loci),
               n_panel_with_loci = sum(panels[[dn]]),
               r_alpha = cor(dd$alpha_outside, dd$alpha_with_loci))
  }))
}))

# The outside fits must be the published ones.
m22 <- fread(MODELS_22)
chk <- merge(tab, m22[, .(donor, model, slope, p)], by = c("donor", "model"))
stopifnot(nrow(chk) == 15)
rob <- fread(ROBUST_22)[test == "alpha ~ pis + bio05 + metro"]
chk2 <- merge(tab[model == "alpha ~ pis + bio05 + metro"],
              rob[, .(donor, slope = estimate, p)], by = "donor")
stopifnot(nrow(chk2) == 3)
chk <- rbind(chk[, .(slope_outside, slope, p_outside, p)],
             chk2[, .(slope_outside, slope, p_outside, p)])
if (!isTRUE(all.equal(chk$slope_outside, chk$slope, tolerance = 1e-8)) ||
    !isTRUE(all.equal(chk$p_outside, chk$p, tolerance = 1e-8)))
  stop("outside-interval fits do not reproduce 22's tables", call. = FALSE)
say("outside-interval fits reproduce %s and the metro row of %s (%d models)",
    basename(MODELS_22), basename(ROBUST_22), nrow(chk))
fwrite(tab, OUT_TABLE)

# ---- console report ----------------------------------------------------
say("\nPanel SNPs and per-pool alpha with and without the three candidate intervals")
for (dn in c("QUI", "PIP", "MOL")) {
  dd <- per_pool[donor == dn]
  say("  %s: %s panel SNPs outside, %s with the intervals; r = %.4f; mean shift %+.4f (range %+.4f to %+.4f)",
      dn, format(sum(panels[[dn]] & !in_loci), big.mark = ","),
      format(sum(panels[[dn]]), big.mark = ","),
      cor(dd$alpha_outside, dd$alpha_with_loci), mean(dd$alpha_diff),
      min(dd$alpha_diff), max(dd$alpha_diff))
}
say("\n%-5s %-28s %-6s %11s %9s %7s   %11s %9s %7s", "donor", "model", "focal",
    "slope out", "P out", "R2 out", "slope with", "P with", "R2 with")
for (i in seq_len(nrow(tab))) {
  r <- tab[i]
  say("%-5s %-28s %-6s %+11.4g %9.3g %7.3f   %+11.4g %9.3g %7.3f", r$donor, r$model,
      r$focal_term, r$slope_outside, r$p_outside, r$adj_r2_outside,
      r$slope_with_loci, r$p_with_loci, r$adj_r2_with_loci)
}

say("\nwrote %s", OUT_TABLE)
say("wrote %s", OUT_PER_POOL)
