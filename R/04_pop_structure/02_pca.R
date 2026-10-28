# =============================================================================
# R/04_pop_structure/02_pca.R
# =============================================================================
# Genome-wide PCA of per-pool allele frequencies across the 27 pools (24 field
# + 3 reference): SNPs with zero variance across pools are dropped, and
# prcomp is run on the scaled pool x SNP matrix. Draws Figure 2 (PC1-PC2, with
# an inset of the field pools labelled by site number) and Figure S2 (PC1-PC3
# and PC2-PC3), and prints the variance explained by PCs 1-5.
#
# Manuscript: Figure 2; Figure S2; Results (PC1 30.2%, PC2 12.2%, PC3 5.9%)
#
# Inputs:  data/sample_categories.csv               (sample metadata)
#          deposit: data/allele_frequencies_v4.csv
# Outputs: figures/final/Figure_2.png
#          figures/final/supplement/Figure_S2.png
#
# Runtime: ~10 min; peak memory several GB (11.6 M SNPs x 27 pools)
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/sample_metadata.R")   # sample_info_all, field_names, ref_names, REFS
source("R/helpers/ref_colors.R")        # ref_cols, ref_labels

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

set.seed(42)   # ggrepel label placement in the Figure 2 inset

AF_FILE <- deposit_file("data/allele_frequencies_v4.csv")
FIG_2   <- "figures/final/Figure_2.png"
FIG_S2  <- "figures/final/supplement/Figure_S2.png"
dir.create(dirname(FIG_S2), recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load allele frequencies + variance filter --------------------------
cat("[", format(Sys.time()), "] Reading allele frequencies:", AF_FILE, "\n")
# Read with fread, selecting only the needed columns. read.csv on this file
# (1.77 GB, 11.6M rows) builds a data.frame several times the file size and,
# with the transpose and prcomp's scaled copy below, exhausts 16 GB of RAM.
# fread + selection gives the same values with far less memory.
keep_cols <- c('chromosome', 'position', 'cumulative_position',
               field_names, ref_names)
avail_cols <- names(data.table::fread(AF_FILE, nrows = 0))
data <- as.data.frame(data.table::fread(AF_FILE,
                                        select = intersect(avail_cols, keep_cols),
                                        showProgress = FALSE))
n_total_snps <- nrow(data)

# Keep metadata columns + the 27 pools (24 field + 3 reference).
data_filtered <- data[, intersect(names(data), keep_cols)]
rm(data); invisible(gc())

# Per-pool AF columns (cols 4 onwards: chromosome, position,
# cumulative_position are the first three metadata cols).
freqs <- data_filtered[, 4:ncol(data_filtered)]

# Drop invariant SNPs (no variance across pools) — these don't contribute
# to PCA and can cause `prcomp(..., scale. = TRUE)` to error on zero-variance
# columns once transposed.
# matrixStats::rowVars gives the same values as apply(freqs, 1, var), computed
# at C level without materialising 11.6M row vectors.
var_per_snp  <- matrixStats::rowVars(as.matrix(freqs), na.rm = TRUE)
n_invariant  <- sum(var_per_snp == 0, na.rm = TRUE)
data_filtered <- data_filtered[var_per_snp > 0, ]

cat("\n[ sanity ]\n")
cat("  Total SNPs in CSV:           ", n_total_snps, "\n")
cat("  SNPs with non-zero variance: ", nrow(data_filtered), "\n")
cat("  Invariant SNPs dropped:      ", n_invariant, "\n")
cat("  Pools retained:              ", ncol(freqs), "\n")
cat("  Pool names: ", paste(colnames(freqs), collapse = ", "), "\n")

# ---- 2. PCA -----------------------------------------------------------------
cat("\n[", format(Sys.time()), "] Running PCA...\n")
rm(freqs); invisible(gc())   # drop the sanity-block copy before prcomp scales

freq_t <- t(as.matrix(data_filtered[, 4:ncol(data_filtered)]))  # samples as rows
rm(data_filtered); invisible(gc())

# rank. = 5 caps the SNP loadings prcomp computes at 5 columns rather than 27
# (the full 11.6M x 27 rotation matrix is ~2.3 GB). Scores and sdev are
# unaffected.
pca_res <- prcomp(freq_t, scale. = TRUE, rank. = 5)
rm(freq_t); invisible(gc())

# Per-PC variance explained (proportion)
variance_explained <- summary(pca_res)$importance[2, ]

# PC scores per pool, joined to metadata
pca_df <- as.data.frame(pca_res$x[, 1:5])
pca_df$sample <- rownames(pca_res$x)
rownames(pca_df) <- NULL
pca_df <- left_join(pca_df, sample_info_all, by = 'sample')

# Site numbers of Figure 1 and Table S1
site_num_df = data.frame('sample' = field_names,
                         'site_num' = c(24,  9,  8, 16,  2,  5, 17,  1, 10,  4, 18, 19,
                                        21, 15,  3, 23, 11, 13,  7, 12, 22,  6, 20, 14))

pca_df = left_join(pca_df, site_num_df, by = 'sample')
fld    <- pca_df[pca_df$sample_type == 'field', ]

# Point fill by pool type: references in their lineage colours, field pools grey.
type_fill <- c(colony_p    = unname(ref_cols['CAL1']),
               field_r     = 'grey50',
               field_s     = 'grey50',
               field_u     = 'grey50',
               reference_p = unname(ref_cols['ME1']),
               reference_q = unname(ref_cols['FL1']))

pc_lab <- function(pc) {
  sprintf("%s (%.1f%% variance)", pc,
          100 * variance_explained[as.integer(sub("PC", "", pc))])
}

# Lineage labels beside each reference point: the italic lineage name over the
# lineage code. side: +1 right of the point, -1 left. vj: vjust of the name
# line (-0.15 name above the point and code below it; -1.45 both above;
# 0.5 name level with the point and code below); the code line sits one line
# lower.
lineage_expr <- c(
  ME1  = 'italic("Cx. pipiens")~"f."~italic("pipiens")',
  FL1  = 'italic("Cx. quinquefasciatus")',
  CAL1 = 'italic("Cx. pipiens")~"f."~italic("molestus")'
)
ref_annotation <- function(df, x, y, side, vj, size, colour = 'black') {
  d <- df[match(REFS, df$sample), c('sample', x, y)]
  names(d) <- c('sample', 'x', 'y')
  s <- unname(side[d$sample])
  d$x       <- d$x + s * 0.025 * diff(range(df[[x]]))
  d$hj      <- ifelse(s > 0, 0, 1)
  d$vj_name <- unname(vj[d$sample])
  d$vj_code <- d$vj_name + 1.25
  d$name    <- unname(lineage_expr[d$sample])
  d$code    <- sprintf('(%s)', ref_labels[d$sample])
  list(
    geom_text(data = d, aes(x = x, y = y, label = name, hjust = hj, vjust = vj_name),
              parse = TRUE, size = size, colour = colour, inherit.aes = FALSE),
    geom_text(data = d, aes(x = x, y = y, label = code, hjust = hj, vjust = vj_code),
              size = size, colour = colour, inherit.aes = FALSE)
  )
}

# ---- 3. Figure 2: PC1-PC2, all pools, with the field-pool inset ------------
# Panel limits are the ggplot default (range + 5% each side) set explicitly, so
# the inset position, given in panel units, converts to data units for the
# dashed lines joining it to the box around the field pools.
expand5 <- function(r) r + c(-1, 1) * 0.05 * diff(r)
xlim2 <- expand5(range(pca_df$PC1))
ylim2 <- expand5(range(pca_df$PC2))

pad   <- 0.015 * diff(range(pca_df$PC1))
box_x <- range(fld$PC1) + c(-pad, pad)
box_y <- range(fld$PC2) + c(-pad, pad)

INSET   <- c(left = 0.48, bottom = 0.10, right = 0.90, top = 0.82)   # panel units
inset_x <- xlim2[1] + unname(INSET[c('left', 'right')]) * diff(xlim2)
inset_y <- ylim2[1] + unname(INSET[c('bottom', 'top')]) * diff(ylim2)

p_fig2 <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  annotate('rect', xmin = box_x[1], xmax = box_x[2], ymin = box_y[1], ymax = box_y[2],
           fill = NA, colour = 'grey20', linewidth = 0.4) +
  annotate('segment', x = box_x[c(2, 1)], y = box_y[c(2, 1)],
           xend = inset_x[c(2, 1)], yend = inset_y[c(2, 1)],
           linetype = 'dashed', colour = 'grey30', linewidth = 0.3) +
  geom_point(aes(fill = sample_type_2),
             shape = 21, color = 'black', size = 6, alpha = 0.4) +
  scale_fill_manual(values = type_fill) +
  annotate('text', x = box_x[2] + pad, y = mean(box_y), label = 'Field\nPopulations',
           hjust = 0, size = 6) +
  ref_annotation(pca_df, 'PC1', 'PC2',
                 side = c(ME1 = -1, FL1 = -1, CAL1 = -1),
                 vj   = c(ME1 = -0.15, FL1 = -0.15, CAL1 = -1.45), size = 6) +
  coord_equal(xlim = xlim2, ylim = ylim2, expand = FALSE) +
  labs(x = pc_lab("PC1"), y = pc_lab("PC2")) +
  theme_bw(base_size = 18) +
  theme(legend.position = 'none')

# Inset: field pools only, labelled by site number
p_fig2_inset <- ggplot(fld, aes(x = PC1, y = PC2)) +
  geom_point(fill = 'grey50', shape = 21, color = 'black', size = 5, alpha = 0.4) +
  geom_text_repel(aes(label = site_num),
                  size = 6, max.overlaps = 10, point.size = 0) +
  coord_equal() +
  theme_bw(base_size = 18) +
  theme(axis.title = element_blank(),
        plot.background = element_rect(fill = 'white', colour = 'black', linewidth = 0.5))

fig2 <- p_fig2 +
  inset_element(p_fig2_inset,
                left  = INSET[['left']],  bottom = INSET[['bottom']],
                right = INSET[['right']], top    = INSET[['top']],
                align_to = 'panel')

ggsave(FIG_2, fig2, width = 12, height = 10, units = 'in', dpi = 300)

# ---- 4. Figure S2: PC1-PC3 and PC2-PC3, all pools --------------------------
s2_panel <- function(x_pc, y_pc, side, vj, field_label_above, tag) {
  gap <- 0.03 * diff(range(pca_df[[y_pc]]))
  if (field_label_above) {   # above the field cluster, ending at its centre
    fl <- list(x = median(fld[[x_pc]]), y = max(fld[[y_pc]]) + gap, hjust = 1, vjust = 0)
  } else {                   # below the field cluster, starting at its left edge
    fl <- list(x = min(fld[[x_pc]]), y = min(fld[[y_pc]]) - gap, hjust = 0, vjust = 1)
  }
  ggplot(pca_df, aes(x = .data[[x_pc]], y = .data[[y_pc]])) +
    geom_point(aes(color = sample_type_2), size = 2) +
    scale_color_manual(values = type_fill) +
    annotate('text', x = fl$x, y = fl$y, label = 'Field\nPopulations',
             hjust = fl$hjust, vjust = fl$vjust, size = 4.2, colour = 'grey30') +
    ref_annotation(pca_df, x_pc, y_pc, side = side, vj = vj,
                   size = 4.2, colour = 'grey30') +
    annotate('text', x = -Inf, y = Inf, label = tag, hjust = -0.6, vjust = 1.4,
             fontface = 'bold', size = 8) +
    coord_equal() +
    labs(x = pc_lab(x_pc), y = pc_lab(y_pc)) +
    theme_bw(base_size = 13) +
    theme(legend.position = 'none')
}

p_s2a <- s2_panel("PC1", "PC3",
                  side = c(ME1 = 1, CAL1 = 1, FL1 = -1),
                  vj   = c(ME1 = 0.5, CAL1 = 0.5, FL1 = -1.45),
                  field_label_above = FALSE, tag = "A")
p_s2b <- s2_panel("PC2", "PC3",
                  side = c(ME1 = -1, CAL1 = 1, FL1 = -1),
                  vj   = c(ME1 = 0.5, CAL1 = 0.5, FL1 = -1.45),
                  field_label_above = TRUE, tag = "B")

ggsave(FIG_S2, p_s2a / p_s2b, width = 8, height = 11.5, units = 'in', dpi = 300)

# ---- 5. Manuscript metrics --------------------------------------------------
cat("\nPCA variance explained (%, PC1-PC5):\n")
print(round(100 * variance_explained[1:5], 2))

cat("\n=== Done. [", format(Sys.time()), "] ===\n")
