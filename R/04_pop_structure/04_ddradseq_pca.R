# =============================================================================
# R/04_pop_structure/04_ddradseq_pca.R
# =============================================================================
# Individual-level PCA of the ddRADseq dataset (A. Arsenault-Benoit and K. L. Bell; pipeline of Bell et al. 2024): posterior
# mean genotypes from entropy version 1 (Gompert et al. 2014; Shastry et al.
# 2021) at 7,258 SNPs, averaged over the two MCMC chains of
# the K = 2 run, centred and not scaled. Pool-seq cannot tell a pool of
# admixed individuals from a mixture of unmixed ones; this ordination answers
# that at individual resolution, and because the FL / ME / CA samples come
# from the source populations of the QUI / PIP / MOL reference pools it also
# places the reference lineages.
#
# Two ordinations are computed:
#   1. all 573 individuals, only to confirm that this pipeline reproduces the
#      original analysis (PC1 5.084%, PC2 2.442%);
#   2. the 553 field and reference individuals, dropping the EV and NF
#      laboratory colonies, which are not from the study region. This is the
#      one plotted.
#
# Manuscript: Figure S4; Supplementary Methods "Individual-level ddRADseq
#             genotyping"
#
# Inputs:  deposit: data/ddradseq/gprob2c0_all.txt
#          deposit: data/ddradseq/gprob2c1_all.txt
#              entropy posterior-mean genotypes, one file per MCMC chain;
#              573 rows (ind_0 .. ind_572) x 7,258 loci plus an id column
#          deposit: data/ddradseq/indiv_ids.txt
#              ind, pop, taxon in the row order of the gprob files. The
#              "taxon" column is the collection-event index (the middle
#              field of the NewID: AS_18_03 -> event 18), not a taxon.
#          deposit: data/ddradseq/ddrad_sample_sheet.csv
#              sample sheet for the 616-individual library: site, PCR species
#              call, urbanisation class, collection date and year. Supplies
#              the wild-caught / colony distinction.
#          Library preparation, alignment, variant calling and the entropy
#          run that produced these files are described in the Supplementary
#          Methods.
# Outputs: figures/final/supplement/Figure_S4.png
#          results/ddradseq_pca_scores.csv   (PC1-PC3 per plotted individual)
#
# Runtime: < 1 min
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/ref_colors.R")   # ref_cols

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

GPROB_C0  <- deposit_file("data/ddradseq/gprob2c0_all.txt")
GPROB_C1  <- deposit_file("data/ddradseq/gprob2c1_all.txt")
IDS_FILE  <- deposit_file("data/ddradseq/indiv_ids.txt")
META_FILE <- deposit_file("data/ddradseq/ddrad_sample_sheet.csv")

FIG_OUT    <- "figures/final/supplement/Figure_S4.png"
SCORES_OUT <- "results/ddradseq_pca_scores.csv"
dir.create(dirname(FIG_OUT), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(SCORES_OUT), showWarnings = FALSE, recursive = TRUE)

# The ddRADseq reference populations map onto our three reference pools:
# FL = Vero Beach (QUI), ME = Orono (PIP), CA = Calumet (MOL). CA is Calumet;
# CL is a field site.
REF_MAP <- c(FL = "QUI", ME = "PIP", CA = "MOL")

# ---------------------------------------------------------------------------
# 1. Load
# ---------------------------------------------------------------------------
message("reading entropy genotype probabilities ...")
c0  <- fread(file = GPROB_C0)
c1  <- fread(file = GPROB_C1)
ids <- fread(file = IDS_FILE)

stopifnot(nrow(c0) == nrow(c1), ncol(c0) == ncol(c1))
stopifnot(nrow(c0) == nrow(ids))
# the two files are chains of one K = 2 run, not two clusters, so they must
# agree on the individual ordering
stopifnot(identical(c0[[1]], c1[[1]]))
stopifnot(identical(c0[[1]], paste0("ind_", seq_len(nrow(c0)) - 1L)))

# EV and NF are laboratory colonies rather than corridor collections. The
# sample sheet is the authority on this, so read it rather than hard-coding.
meta <- fread(file = META_FILE)
ids <- merge(ids, meta[, .(NewID, Location, ID, Prim_Class, Year)],
             by.x = "ind", by.y = "NewID", all.x = TRUE, sort = FALSE)
setkey(ids, NULL)
stopifnot(!anyNA(ids$Location))

# ---------------------------------------------------------------------------
# 2. Average the chains
# ---------------------------------------------------------------------------
G <- (as.matrix(c0[, -1]) + as.matrix(c1[, -1])) / 2
rm(c0, c1); invisible(gc())
message(sprintf("genotype matrix: %d individuals x %d loci", nrow(G), ncol(G)))

pct <- function(p) 100 * (p$sdev^2 / sum(p$sdev^2))

run_pca <- function(mat) {
  keep <- apply(mat, 2, function(x) stats::var(x) > 0)
  if (any(!keep)) message(sprintf("  dropping %d zero-variance loci", sum(!keep)))
  # centred, not scaled: this is what reproduces the original analysis
  prcomp(mat[, keep, drop = FALSE], scale. = FALSE, rank. = 5)
}

# ---------------------------------------------------------------------------
# 3. Provenance check on the full set
# ---------------------------------------------------------------------------
message("full set (573) -- check against the original analysis:")
pca_all <- run_pca(G)
ve_all <- pct(pca_all)
message(sprintf("  PC1 %.3f%%  PC2 %.3f%%  (original: 5.084%%, 2.442%%)",
                ve_all[1], ve_all[2]))
if (max(abs(ve_all[1:2] - c(5.084, 2.442))) > 0.01) {
  warning("full-set percentages no longer match the original analysis")
}

# ---------------------------------------------------------------------------
# 4. The plotted ordination: drop the colony populations
# ---------------------------------------------------------------------------
is_colony <- tolower(ids$Prim_Class) == "colony" & !(ids$pop %in% names(REF_MAP))
message(sprintf("dropping %d colony individuals (%s)", sum(is_colony),
                paste(sort(unique(ids$pop[is_colony])), collapse = ", ")))

ids_sub <- ids[!is_colony]
pca_res <- run_pca(G[!is_colony, , drop = FALSE])
ve <- pct(pca_res)
message(sprintf("plotted set (%d) -- PC1 %.3f%%  PC2 %.3f%%",
                nrow(ids_sub), ve[1], ve[2]))

# The sign of a principal component is arbitrary. Orient PC1 so the
# quinquefasciatus reference sits on the right, matching Figure 2.
if (mean(pca_res$x[ids_sub$pop == "FL", 1]) < 0) {
  pca_res$x[, 1] <- -pca_res$x[, 1]
  message("flipped PC1 so QUI is positive, matching Figure 2")
}

# ---------------------------------------------------------------------------
# 5. Assemble the plotting frame
# ---------------------------------------------------------------------------
df <- data.table(
  ind = ids_sub$ind,
  pop = ids_sub$pop,
  site = ids_sub$Location,
  pcr_id = ids_sub$ID,
  year = ids_sub$Year,
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  PC3 = pca_res$x[, 3]
)
df[, group := fifelse(pop %in% names(REF_MAP), unname(REF_MAP[pop]), "Field")]
df[, group := factor(group, levels = c("Field", "QUI", "PIP", "MOL"))]

message("group sizes: ", paste(sprintf("%s=%d", levels(df$group), table(df$group)),
                               collapse = "  "))

# Check of the claim the figure is used to make: the quinquefasciatus
# reference should sit apart, with empty space between it and everything else.
qui_rng   <- range(df[group == "QUI", PC1])
other_rng <- range(df[group != "QUI", PC1])
message(sprintf("PC1: QUI [%.2f, %.2f] vs all others [%.2f, %.2f] -- empty gap = %.2f",
                qui_rng[1], qui_rng[2], other_rng[1], other_rng[2],
                min(qui_rng) - max(other_rng)))

fwrite(df, SCORES_OUT)

# ---------------------------------------------------------------------------
# 6. Figure S4
# ---------------------------------------------------------------------------
grp_cols <- c(
  Field = "grey65",
  QUI   = unname(ref_cols["FL1"]),
  PIP   = unname(ref_cols["ME1"]),
  MOL   = unname(ref_cols["CAL1"])
)

p <- ggplot(df, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = group), shape = 21, colour = "black",
             size = 2.6, stroke = 0.25, alpha = 0.9) +
  scale_fill_manual(values = grp_cols, name = NULL,
                    breaks = c("Field", "QUI", "PIP", "MOL"),
                    labels = c("Field populations", "QUI", "PIP", "MOL")) +
  labs(x = sprintf("PC1 (%.1f%% variance)", ve[1]),
       y = sprintf("PC2 (%.1f%% variance)", ve[2])) +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave(FIG_OUT, p, width = 9, height = 7, dpi = 300)

message("-> ", FIG_OUT)
