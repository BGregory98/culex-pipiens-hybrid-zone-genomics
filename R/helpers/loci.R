# =============================================================================
# R/helpers/loci.R
# =============================================================================
#
# Coordinates of the candidate loci on idCulPipi1 (1-based, inclusive). Every
# script that needs a locus interval sources this file instead of typing the
# numbers, so the paper has one definition of each locus.
#
#   source("R/helpers/loci.R")   # working directory = repository root
#
# CpL3a is the chromosome 3 inversion between the DELLY breakpoints in the MOL
# reference pool (CAL1). DELLY reports the two junction orientations of the
# inversion as separate high-confidence calls:
#   3to3  149,402,130 - 161,021,387  (paired-end support 6, MAPQ 44)
#   5to5  149,402,438 - 161,021,767  (paired-end support 10, MAPQ 60)
# CpL3a is the interval bounded by both: 149,402,438 - 161,021,387 (11.62 Mb).
# The genome-wide null (cluster/07_multilocus/01_genome_wide_null.R) uses the
# same interval; its chr3 random regions are this width.
# =============================================================================

CHR1 <- "ENA_OZ004311_OZ004311.1"
CHR2 <- "ENA_OZ004312_OZ004312.1"
CHR3 <- "ENA_OZ004313_OZ004313.1"

CPL1A_START <- 108500000L; CPL1A_END <- 111000000L   # CYP6 haplotype block
CPL2A_START <-  12500000L; CPL2A_END <-  15000000L   # Yan-centred haplotype block
CPL3A_START <- 149402438L; CPL3A_END <- 161021387L   # inversion (DELLY, CAL1)

# The 500 kb core of CpL2a around Yan, where the BayPass candidates sit. The
# clines, f3 windows, LD sets, Table S4 panels and the panel-based tests use the
# core; the genome-wide null, the per-pool ancestry test and the haplotype-block
# diagnostics use the 2.5 Mb block above.
CPL2A_CORE_START <- 13500000L; CPL2A_CORE_END <- 14000000L

LOCI <- list(
  CpL1a = list(chr = CHR1, start = CPL1A_START, end = CPL1A_END),
  CpL2a = list(chr = CHR2, start = CPL2A_START, end = CPL2A_END),
  CpL3a = list(chr = CHR3, start = CPL3A_START, end = CPL3A_END)
)

# Genes named in the paper (lifted Cx. quinquefasciatus JHB annotation IDs).
# The six CYP6 paralogs of the CpL1a cluster, in chromosome order.
CYP6_GENES <- c("CQUJHB009005", "CQUJHB008986", "CQUJHB009148",
                "CQUJHB017007", "CQUJHB001163", "CQUJHB007871")
CYP6_NAMES <- c(CQUJHB009005 = "CYP6AA7-like", CQUJHB008986 = "CYP6AA8-like",
                CQUJHB009148 = "CYP6AA9-like", CQUJHB017007 = "CYP6BZ2-like",
                CQUJHB001163 = "CYP6CC2-like", CQUJHB007871 = "CYP6BB4-like")   # P450Atlas
YAN_GENE   <- "CQUJHB014704"   # CpL2a
VGSC_GENE  <- "CQUJHB005371"   # voltage-gated sodium channel (para), kdr site
