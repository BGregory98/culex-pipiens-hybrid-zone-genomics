# =============================================================================
# R/07_multi_locus/21_private_panel_sizes.R
# =============================================================================
# Donor-private panel sizes in each of the three candidate intervals and in the
# remaining genomic background: the full 4 x 3 matrix of which comparisons the
# data can and cannot support. Panels are defined as everywhere else: an allele
# at frequency >= 0.95 in one reference pool and <= 0.05 in both others, counted
# among positions with a called frequency in all three references.
#
# Manuscript: Table S4; Results "Corridor populations are weakly differentiated
#             and predominantly Cx. p. f. pipiens-like" (genome-wide panel sizes)
#             and "Measuring donor ancestry and testing it against a genome-wide
#             null"; Supplementary Results "Background composition estimated from
#             lineage-private panels", "Reference-pool provenance and its effect
#             on the private panels", "Donor-private panels give opposing ancestry
#             gradients at CpL1a and CpL2a" and "Locus statistics against a
#             differentiation-matched null" (panel enrichment)
#
# Inputs:  deposit: data/allele_frequencies_v4.csv   (R/01_filtering/01_allele_frequency_table.R)
# Outputs: results/multilocus/Table_panel_sizes_by_region.csv
#
# Runtime: seconds
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")
suppressMessages(library(data.table))

AF_PATH <- deposit_file("data/allele_frequencies_v4.csv")

HI <- 0.95; LO <- 0.05
# CpL2a is the 500 kb candidate core (13.5-14.0 Mb) in the private-panel analyses.
LOCI <- list(
  CpL1a = list(chr = CHR1, start = CPL1A_START, end = CPL1A_END),
  CpL2a = list(chr = CHR2, start = CPL2A_CORE_START,      end = CPL2A_CORE_END),
  CpL3a = list(chr = CHR3, start = CPL3A_START, end = CPL3A_END)
)

say("reading %s ...", AF_PATH)
af <- fread(AF_PATH, select = c("chromosome", "position", REFS), showProgress = FALSE)
say("matrix rows: %s", format(nrow(af), big.mark = ","))

q <- af[[QUI]]; p <- af[[PIP]]; m <- af[[MOL]]
ok <- is.finite(q) & is.finite(p) & is.finite(m)
PANEL <- list(
  QUI = ok & q >= HI & p <= LO & m <= LO,
  PIP = ok & p >= HI & q <= LO & m <= LO,
  MOL = ok & m >= HI & q <= LO & p <= LO
)

in_cand <- rep(FALSE, nrow(af))
region <- rep("Background", nrow(af))
for (nm in names(LOCI)) {
  L <- LOCI[[nm]]
  k <- af$chromosome == L$chr & af$position >= L$start & af$position <= L$end
  region[k] <- nm
  in_cand <- in_cand | k
}

rows <- list()
for (rg in c("CpL1a", "CpL2a", "CpL3a", "Background")) {
  k <- region == rg
  rows[[rg]] <- data.table(
    region      = rg,
    total_SNPs  = sum(k & ok),
    QUI_private = sum(PANEL$QUI & k),
    PIP_private = sum(PANEL$PIP & k),
    MOL_private = sum(PANEL$MOL & k))
}
tab <- rbindlist(rows)

say("")
say("=== donor-private panel sizes, by region (>= %.2f / <= %.2f) ===", HI, LO)
say("%-11s %12s %12s %12s %12s", "region", "SNPs", "QUI-private", "PIP-private", "MOL-private")
for (i in seq_len(nrow(tab))) {
  say("%-11s %12s %12s %12s %12s", tab$region[i],
      format(tab$total_SNPs[i], big.mark = ","),
      format(tab$QUI_private[i], big.mark = ","),
      format(tab$PIP_private[i], big.mark = ","),
      format(tab$MOL_private[i], big.mark = ","))
}
say("")
say("genome-wide totals: QUI %s, PIP %s, MOL %s",
    format(sum(PANEL$QUI), big.mark = ","),
    format(sum(PANEL$PIP), big.mark = ","),
    format(sum(PANEL$MOL), big.mark = ","))
say("")
say("Reading: PIP is the least represented panel everywhere, which is expected -- a")
say("reference with the highest heterozygosity is least often fixed for an allele")
say("absent from both others. At CpL1a it falls to a single SNP, which is why no")
say("PIP ancestry estimate is reported there.")

fwrite(tab, "results/multilocus/Table_panel_sizes_by_region.csv")
say("")
say("wrote results/multilocus/Table_panel_sizes_by_region.csv")
