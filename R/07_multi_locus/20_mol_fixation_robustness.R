# =============================================================================
# R/07_multi_locus/20_mol_fixation_robustness.R
# =============================================================================
# Two robustness checks.
#   1. Metropolitan assignment. Each site is assigned to the nearer city centre
#      by great-circle distance. Lists the sites within 10 km of equal distance,
#      recomputes the within-group bio05-impervious surface correlation under
#      alternative rules (latitude splits; ambiguous sites reassigned), and gives
#      the bio05 range carried within each metropolitan area.
#   2. MOL fixation at CpL3a. The MOL-private panel is the largest of the three
#      and carries CpL3a. A "MOL-private" SNP where the corridor mean is already
#      high may reflect MOL fixation on a widespread allele rather than
#      lineage-restricted variation, so the corridor-wide median Delta is
#      recomputed while excluding panel SNPs above a corridor-mean cutoff.
#
# Manuscript: Supplementary Results "Temperature is collinear with position along
#             the corridor" (metropolitan-assignment paragraph; within-area bio05
#             spans); Supplementary Results "Robustness of the CpL3a MOL-like
#             signal to MOL fixation"; Results CpL3a section ("robust to the
#             exclusion of corridor-common alleles")
#
# Inputs:  data/site_coords_with_spatial_data_allvars.csv   (R/00_setup/01_site_predictors.R)
#          deposit: data/allele_frequencies_v4.csv          (R/01_filtering/01_allele_frequency_table.R)
# Outputs: none (console only)
#
# Runtime: under a minute
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")
source("R/helpers/geo.R")
suppressMessages(library(data.table))

AF_PATH <- deposit_file("data/allele_frequencies_v4.csv")
need_files("data/site_coords_with_spatial_data_allvars.csv")

CPL3A <- list(chr = CHR3, start = CPL3A_START, end = CPL3A_END)
# CpL2a is the 500 kb candidate core (13.5-14.0 Mb) in the private-panel analyses.
LOCI <- list(
  CpL1a = list(chr = CHR1, start = CPL1A_START, end = CPL1A_END),
  CpL2a = list(chr = CHR2, start = CPL2A_CORE_START,      end = CPL2A_CORE_END),
  CpL3a = CPL3A
)

# -----------------------------------------------------------------------------
# Sites and metropolitan assignment
# -----------------------------------------------------------------------------
env <- fread("data/site_coords_with_spatial_data_allvars.csv", showProgress = FALSE)
setkey(env, site); e <- env[field_names]

e$d_bal <- hav(e$lat, e$lon, BAL[1], BAL[2])
e$d_dc  <- hav(e$lat, e$lon, DC[1],  DC[2])
e$metro <- ifelse(e$d_bal < e$d_dc, "Baltimore", "Washington")
e$margin <- abs(e$d_bal - e$d_dc)

# -----------------------------------------------------------------------------
# 1a. How much do the within-metro correlations depend on the rule?
# -----------------------------------------------------------------------------
say("=== the assignment rule, and its ambiguous cases ===")
say("Rule: each site is assigned to whichever city centre is nearer by great-circle")
say("distance -- Baltimore (39.2904 N, 76.6122 W) or Washington (38.9072 N, 77.0369 W).")
say("split: Baltimore %d, Washington %d", sum(e$metro == "Baltimore"),
    sum(e$metro == "Washington"))
say("")
say("sites within 10 km of equal distance (the ambiguous mid-corridor set):")
amb <- e[margin < 10][order(margin)]
if (nrow(amb)) print(amb[, .(site, metro, d_bal = round(d_bal, 1), d_dc = round(d_dc, 1),
                             margin_km = round(margin, 2), bio05 = round(bio05, 3))],
                     row.names = FALSE) else say("  none")

wm <- function(grp, lab) {
  r <- sapply(unique(grp), function(g) {
    k <- grp == g
    if (sum(k) < 4) return(NA_real_)
    cor(e$bio05[k], e$pis[k])
  })
  n <- sapply(unique(grp), function(g) sum(grp == g))
  say("  %-34s %s", lab,
      paste(sprintf("%s n=%d r=%+.3f", unique(grp), n, r), collapse = " | "))
}
say("")
say("r(bio05, pis) within group, under alternative rules:")
wm(e$metro, "nearest city centre (used)")
wm(ifelse(e$lat >= 39.10, "North", "South"), "latitude split at 39.10 N")
wm(ifelse(e$lat >= 39.00, "North", "South"), "latitude split at 39.00 N")
wm(ifelse(e$lat >= 39.20, "North", "South"), "latitude split at 39.20 N")
# nearest-centre, but forcing every ambiguous site the other way
alt <- e$metro; alt[e$margin < 10] <- ifelse(alt[e$margin < 10] == "Baltimore",
                                            "Washington", "Baltimore")
wm(alt, "nearest centre, ambiguous flipped")
say("")
say("Reading: latitude splits are NOT the rule and do not reproduce the Baltimore")
say("value; the corridor runs southwest-northeast so distance and latitude differ.")

# -----------------------------------------------------------------------------
# 1b. Thermal range carried within each metropolitan area
# -----------------------------------------------------------------------------
say("")
say("=== thermal range carried within each metropolitan area ===")
say("corridor-wide bio05 range: %.3f to %.3f C (span %.3f)",
    min(e$bio05), max(e$bio05), diff(range(e$bio05)))
for (mm in c("Washington", "Baltimore")) {
  k <- e$metro == mm
  say("  %-11s n=%2d  bio05 %.3f to %.3f C (span %.3f = %.0f%% of corridor)",
      mm, sum(k), min(e$bio05[k]), max(e$bio05[k]), diff(range(e$bio05[k])),
      100 * diff(range(e$bio05[k])) / diff(range(e$bio05)))
}

# -----------------------------------------------------------------------------
# 2. Does CpL3a survive dropping MOL-fixed-but-corridor-high SNPs?
# -----------------------------------------------------------------------------
hdr   <- names(fread(AF_PATH, nrows = 0, showProgress = FALSE))
pools <- setdiff(hdr, c("chromosome", "position", "cumulative_position"))
fieldp <- setdiff(pools, REFS)
say("")
say("reading %s ...", AF_PATH)
af <- fread(AF_PATH, select = c("chromosome", "position", pools), showProgress = FALSE)
q <- af[[QUI]]; p <- af[[PIP]]; m <- af[[MOL]]
ok <- is.finite(q) & is.finite(p) & is.finite(m)
in_cand <- rep(FALSE, nrow(af))
for (L in LOCI) in_cand <- in_cand | (af$chromosome == L$chr &
                                      af$position >= L$start & af$position <= L$end)
F <- as.matrix(af[, ..fieldp])
corridor <- rowMeans(F, na.rm = TRUE)

mol_panel <- ok & m >= 0.95 & q <= 0.05 & p <= 0.05
loc <- mol_panel & af$chromosome == CPL3A$chr &
       af$position >= CPL3A$start & af$position <= CPL3A$end
gen <- mol_panel & !in_cand

say("")
say("=== CpL3a MOL-private Delta, excluding MOL-fixed / corridor-high SNPs ===")
say("A 'MOL-private' SNP where the corridor mean is already high may reflect MOL")
say("fixation on a widespread allele rather than lineage-restricted variation.")
say("%-28s %8s %8s %10s %10s", "corridor-mean exclusion", "n_locus", "n_gen",
    "median a", "median D")
for (cut in c(1.00, 0.80, 0.60, 0.50, 0.40, 0.30, 0.20)) {
  keep <- corridor <= cut
  l <- loc & keep; g <- gen & keep
  if (sum(l) < 20 || sum(g) < 100) {
    say("%-28s %8d %8d   too few SNPs retained",
        sprintf("drop corridor mean > %.2f", cut), sum(l), sum(g)); next
  }
  a_l <- colMeans(F[l, , drop = FALSE], na.rm = TRUE)
  a_g <- colMeans(F[g, , drop = FALSE], na.rm = TRUE)
  say("%-28s %8s %8s %10.4f %+10.4f",
      if (cut >= 1) "none (all panel SNPs)" else sprintf("drop corridor mean > %.2f", cut),
      format(sum(l), big.mark = ","), format(sum(g), big.mark = ","),
      median(a_l), median(a_l - a_g))
}
say("")
say("Reading: if median Delta holds as the exclusion tightens, the CpL3a signal is")
say("not an artifact of MOL being fixed for alleles that are common corridor-wide.")
