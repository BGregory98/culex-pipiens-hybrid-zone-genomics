# =============================================================================
# R/08_loci/06_cyp6_microsynteny.R
# =============================================================================
# Draws the microsynteny of the CpL1a CYP6 cluster between chromosome 1 of the
# Cx. quinquefasciatus JHB2020 assembly and the Cx. pipiens idCulPipi1
# reference, from a minimap2 (-x asm10) alignment of the two chromosomes. The
# two assemblies align in opposite orientation, so JHB2020 coordinates are
# reversed. A window of each assembly spanning the six paralogs (+/- 5 kb) is
# shown as (A) a gene-order schematic and (B) a dot-plot of the primary
# alignments (length >= 1 kb, MAPQ >= 5) between the two windows.
#
# Manuscript: Figure S9; Supplementary Results "The CpL1a CYP6 cluster:
#             phylogenetic placement, microsynteny and sequence similarity";
#             Supplementary Methods "Coding consequences, phylogenetic
#             placement, microsynteny and consensus sequences of the CpL1a CYP6
#             paralogs"
#
# Inputs:  deposit: data/reference_genome/chr1_only/chr1_quinq_vs_pip_chr1.asm10.paf   (cluster/08_phylogeny/chr1_minimap.sh)
#          deposit: data/reference_genome/chr1_only/chr1_quinq.gff3
#          deposit: data/reference_genome/chr1_only/chr1_pipiens.gff3
# Outputs: figures/final/supplement/Figure_S9.png
#
# Runtime: < 1 min
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(tibble)
library(patchwork)

# ---- files ----
paf_file    <- deposit_file("data/reference_genome/chr1_only/chr1_quinq_vs_pip_chr1.asm10.paf")
quinq_gff   <- deposit_file("data/reference_genome/chr1_only/chr1_quinq.gff3")
pipiens_gff <- deposit_file("data/reference_genome/chr1_only/chr1_pipiens.gff3")

FIG_OUT <- "figures/final/supplement/Figure_S9.png"
dir.create(dirname(FIG_OUT), recursive = TRUE, showWarnings = FALSE)

# ---- whole-chromosome alignment filters ----
min_aln_len <- 100
min_mapq <- 5

# ---- how much extra sequence to show around the focal cluster ----
padding_bp <- 5000

# ---- helper: extract a GFF attribute ----
get_attr <- function(x, key) {
  stringr::str_match(x, paste0("(?:^|;)", key, "=([^;]+)"))[, 2]
}

# ---- read PAF ----
paf_raw <- readr::read_delim(
  paf_file,
  delim = "\t",
  col_names = FALSE,
  quote = "",
  show_col_types = FALSE
)

if (ncol(paf_raw) < 12) {
  stop("PAF appears to have fewer than 12 required columns.")
}

if (ncol(paf_raw) > 12) {
  tag_text <- apply(as.data.frame(paf_raw[, -(1:12), drop = FALSE]), 1, paste, collapse = "\t")
} else {
  tag_text <- rep("", nrow(paf_raw))
}

paf <- paf_raw %>%
  transmute(
    qname  = X1,  # quinq
    qlen   = as.numeric(X2),
    qstart = as.numeric(X3) + 1,
    qend   = as.numeric(X4),
    strand = X5,
    tname  = X6,  # pipiens
    tlen   = as.numeric(X7),
    tstart = as.numeric(X8) + 1,
    tend   = as.numeric(X9),
    nmatch = as.numeric(X10),
    alnlen = as.numeric(X11),
    mapq   = as.numeric(X12),
    tags   = tag_text
  )

# keep primary alignments if tags are present
if (any(str_detect(paf$tags, "tp:A:P"))) {
  paf <- paf %>%
    mutate(primary = str_detect(tags, "tp:A:P")) %>%
    filter(primary)
}

paf <- paf %>%
  filter(alnlen >= min_aln_len, mapq >= min_mapq)

if (nrow(paf) == 0) {
  stop("No PAF alignments remain after filtering. Relax min_aln_len or min_mapq.")
}

if (length(unique(paf$qname)) != 1 || length(unique(paf$tname)) != 1) {
  stop("Expected chr1-only PAF with one query seq and one target seq.")
}

q_len <- unique(paf$qlen)
t_len <- unique(paf$tlen)

if (length(q_len) != 1 || length(t_len) != 1) {
  stop("Expected exactly one chromosome length per genome in chr1-only PAF.")
}

q_len <- q_len[1]
t_len <- t_len[1]

# flip quinq so local order matches pipiens left-to-right
paf <- paf %>%
  mutate(
    qstart_plot = q_len - qend + 1,
    qend_plot   = q_len - qstart + 1,
    pip_mid     = (tstart + tend) / 2,
    quinq_mid   = (qstart_plot + qend_plot) / 2
  )

# ---- read GFF ----
read_gff_genes <- function(file, seqid_keep, genome_label, y_pos) {
  readr::read_delim(
    file,
    delim = "\t",
    comment = "#",
    quote = "",
    col_names = c("seqid","source","type","start","end","score","strand","phase","attributes"),
    col_types = cols(
      seqid      = col_character(),
      source     = col_character(),
      type       = col_character(),
      start      = col_double(),
      end        = col_double(),
      score      = col_character(),
      strand     = col_character(),
      phase      = col_character(),
      attributes = col_character()
    )
  ) %>%
    filter(seqid == seqid_keep, type %in% c("gene", "protein_coding_gene")) %>%
    mutate(
      genome = genome_label,
      y = y_pos,
      gene_id = get_attr(attributes, "ID"),
      gene_name = coalesce(
        get_attr(attributes, "Name"),
        get_attr(attributes, "gene"),
        get_attr(attributes, "gene_name"),
        gene_id
      ),
      label = gene_name
    )
}

genes_quinq <- read_gff_genes(quinq_gff, "chr1_quinq", "quinq", 2)
genes_pip   <- read_gff_genes(pipiens_gff, "chr1_pipiens", "pipiens", 1)

if (nrow(genes_quinq) == 0) stop("No quinq genes loaded.")
if (nrow(genes_pip) == 0) stop("No pipiens genes loaded.")

genes_quinq <- genes_quinq %>%
  mutate(
    start_plot = q_len - end + 1,
    end_plot   = q_len - start + 1,
    focal      = gene_id %in% CYP6_GENES | gene_name %in% CYP6_GENES
  )

genes_pip <- genes_pip %>%
  mutate(
    start_plot = start,
    end_plot   = end,
    focal      = gene_id %in% CYP6_GENES | gene_name %in% CYP6_GENES
  )

# ---- define focal windows directly from focal genes in both genomes ----
focal_pip <- genes_pip %>% filter(focal)
focal_quinq <- genes_quinq %>% filter(focal)

if (nrow(focal_pip) == 0) {
  stop("No focal genes found on pipiens track. Check the CYP6 gene IDs.")
}
if (nrow(focal_quinq) == 0) {
  stop("No focal genes found on quinq track. Check the CYP6 gene IDs.")
}

pipiens_window <- c(
  max(1, min(focal_pip$start_plot) - padding_bp),
  min(t_len, max(focal_pip$end_plot) + padding_bp)
)

quinq_window <- c(
  max(1, min(focal_quinq$start_plot) - padding_bp),
  min(q_len, max(focal_quinq$end_plot) + padding_bp)
)

# ---- subset genes to local windows ----
genes_pip_plot <- genes_pip %>%
  filter(end_plot >= pipiens_window[1], start_plot <= pipiens_window[2]) %>%
  mutate(
    start_local = start_plot - pipiens_window[1],
    end_local   = end_plot   - pipiens_window[1]
  )

genes_quinq_plot <- genes_quinq %>%
  filter(end_plot >= quinq_window[1], start_plot <= quinq_window[2]) %>%
  mutate(
    start_local = start_plot - quinq_window[1],
    end_local   = end_plot   - quinq_window[1]
  )

genes_plot <- bind_rows(genes_pip_plot, genes_quinq_plot)

# P450Atlas top matches (the "-like" labels used in the text)
genes_plot$label[genes_plot$label == 'CQUJHB009005'] = 'CYP6AA7'
genes_plot$label[genes_plot$label == 'CQUJHB008986'] = 'CYP6AA8'
genes_plot$label[genes_plot$label == 'CQUJHB009148'] = 'CYP6AA9'
genes_plot$label[genes_plot$label == 'CQUJHB017007'] = 'CYP6BZ2'
genes_plot$label[genes_plot$label == 'CQUJHB001163'] = 'CYP6CC2'
genes_plot$label[genes_plot$label == 'CQUJHB007871'] = 'CYP6BB4'

# ---- local alignments for the window view ----
# A tens-of-kb window needs a more permissive filter than the whole-chromosome
# view: primary hits of length >= 1 kb and MAPQ >= 5.
min_local_aln_len <- 1000
min_local_mapq    <- 5

paf_local <- paf_raw %>%
  transmute(
    qname  = X1,  # quinq
    qlen   = as.numeric(X2),
    qstart = as.numeric(X3) + 1,
    qend   = as.numeric(X4),
    strand = X5,
    tname  = X6,  # pipiens
    tlen   = as.numeric(X7),
    tstart = as.numeric(X8) + 1,
    tend   = as.numeric(X9),
    nmatch = as.numeric(X10),
    alnlen = as.numeric(X11),
    mapq   = as.numeric(X12),
    tags   = if (ncol(paf_raw) > 12) {
      apply(as.data.frame(paf_raw[, -(1:12), drop = FALSE]), 1, paste, collapse = "\t")
    } else {
      rep("", nrow(paf_raw))
    }
  )

# keep primary alignments if tp:A:P is present
if (any(str_detect(paf_local$tags, "tp:A:P"))) {
  paf_local <- paf_local %>%
    mutate(primary = str_detect(tags, "tp:A:P")) %>%
    filter(primary)
}

# local plotting filter
paf_local <- paf_local %>%
  filter(alnlen >= min_local_aln_len, mapq >= min_local_mapq) %>%
  mutate(
    # flip quinq so both genomes are plotted in the same left-to-right orientation
    qstart_plot = q_len - qend + 1,
    qend_plot   = q_len - qstart + 1,

    # local coordinates inside the focal windows
    pip_start_local   = tstart       - pipiens_window[1],
    pip_end_local     = tend         - pipiens_window[1],
    quinq_start_local = qstart_plot  - quinq_window[1],
    quinq_end_local   = qend_plot    - quinq_window[1],

    pip_mid_local     = (tstart + tend) / 2      - pipiens_window[1],
    quinq_mid_local   = (qstart_plot + qend_plot) / 2 - quinq_window[1]
  ) %>%
  filter(
    tend >= pipiens_window[1], tstart <= pipiens_window[2],
    qend_plot >= quinq_window[1], qstart_plot <= quinq_window[2]
  )

if (nrow(paf_local) == 0) {
  stop("No local alignments remain. Lower min_local_aln_len/min_local_mapq or re-align the local windows.")
}

message("Local PAF rows plotted: ", nrow(paf_local))

# ---- local window sizes ----
pip_local_len   <- pipiens_window[2] - pipiens_window[1]
quinq_local_len <- quinq_window[2]   - quinq_window[1]
local_xmax      <- max(pip_local_len, quinq_local_len)

# ---- (A) gene-order schematic ----
chr_bars <- tibble(
  genome = c("pipiens", "quinq"),
  y      = c(1, 2),
  label  = c('italic("Cx. pipiens")~"assembly (idCulPipi1)"',
             'italic("Cx. quinquefasciatus")~"assembly (JHB2020)"')
)

p_order <- ggplot() +
  geom_rect(
    data = chr_bars,
    aes(xmin = -Inf, xmax = Inf, ymin = y - 0.06, ymax = y + 0.06),
    fill = "grey70",
    color = "grey30",
    linewidth = 0.2
  ) +
  geom_rect(
    data = genes_plot %>% filter(!focal),
    aes(
      xmin = start_local,
      xmax = end_local,
      ymin = y - 0.06,
      ymax = y + 0.06
    ),
    fill = "grey70",
    color = "black",
    linewidth = 0.2
  ) +
  geom_rect(
    data = genes_plot %>% filter(focal),
    aes(
      xmin = start_local,
      xmax = end_local,
      ymin = y - 0.08,
      ymax = y + 0.08
    ),
    fill = "firebrick",
    color = "black",
    linewidth = 0.3
  ) +
  geom_text(
    data = genes_plot %>% filter(focal),
    aes(
      x = (start_local + end_local) / 2,
      y = y + 0.1,
      label = label
    ),
    angle = 45,
    hjust = 0,
    size = 3
  ) +
  geom_text(
    data = chr_bars,
    aes(x = -Inf, y = y - 0.17, label = label),
    parse = TRUE,
    hjust = -0.03,
    size = 3.9
  ) +
  annotate("text", x = -Inf, y = Inf, label = "A", fontface = "bold",
           size = 6.5, hjust = -0.6, vjust = 1.4) +
  scale_y_continuous(limits = c(0.72, 2.28)) +
  coord_cartesian(xlim = c(0, local_xmax)) +
  labs(x = "Local position within syntenic window (bp)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---- (B) local dot-plot ----
# x = pipiens local coordinate
# y = quinq local coordinate (flipped so conserved synteny should trend positive)
p_dot <- ggplot(paf_local) +
  geom_segment(
    aes(
      x = pip_start_local,
      xend = pip_end_local,
      y = quinq_start_local,
      yend = quinq_end_local
    ),
    linewidth = 0.5,
    alpha = 0.8
  ) +
  annotate("text", x = -Inf, y = Inf, label = "B", fontface = "bold",
           size = 6.5, hjust = -0.6, vjust = 1.4) +
  coord_equal(
    xlim = c(0, pip_local_len),
    ylim = c(0, quinq_local_len),
    expand = FALSE
  ) +
  labs(
    x = expression(italic("Cx. pipiens")~"local position (bp)"),
    y = expression(atop(italic("Cx. quinquefasciatus")~"local position", "(bp, flipped)"))
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank()
  )

fig <- p_order / p_dot

ggsave(FIG_OUT, fig, width = 7, height = 7, units = "in", dpi = 300)
message("wrote ", FIG_OUT)
