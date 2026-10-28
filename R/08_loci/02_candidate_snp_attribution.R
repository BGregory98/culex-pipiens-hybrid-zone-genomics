# =============================================================================
# R/08_loci/02_candidate_snp_attribution.R
# =============================================================================
# Attributes the BayPass candidate SNPs to the protein-coding genes of the
# lifted JHB2020 annotation in three intervals: CpL1a (every gene, so the six
# CYP6 paralogs are ranked against the rest of the interval), the 500 kb CpL2a
# core around Yan (13.5-14.0 Mb), and the CpL3a inversion. For each gene it
# counts candidate SNPs in the gene span +/- 2 kb and in the isoform-union CDS,
# with densities, the maximum XtXst and the fraction on the locus's nominal
# predictor. Also lists the candidate SNPs in and around Vgsc (para).
#
# Manuscript: Tables S5, S7, S8; Results "Cx. quinquefasciatus ancestry and
#             sequencing read depth at a CYP6 cluster rise with summer
#             temperature" (79 genes, CYP6 ranks), "Cx. quinquefasciatus ancestry
#             is excluded at a developmental-quiescence locus on the same
#             gradient" (Yan, 17 of 110); Supplementary Results "Candidate-SNP
#             attribution within the CpL1a interval", "The CpL3a inversion
#             shows no localized gene-level candidate-SNP concentration";
#             Supplementary Methods "Annotation liftover, gene-content surveys
#             and candidate-SNP attribution"
#
# Inputs:  results/baypass/v4/baypass_outliers_clustered.csv   (R/03_baypass/01_baypass_outliers.R)
#          deposit: data/quinq_on_pip.gff3_polished
# Outputs: results/chr1/Table_chr1_cds_snp_attribution.csv     (Table S5)
#          results/chr2/Table_chr2_cds_snp_attribution.csv     (Table S7)
#          results/chr3/Table_chr3_cds_snp_attribution.csv     (Table S8)
#
# Runtime: ~2 min (GFF import)
# =============================================================================
#
# Why multi-metric, not just CDS density. Yan's candidate SNPs are almost all
# intronic or UTR. Adaptive variation in insects frequently targets
# cis-regulatory and intronic regions rather than coding sequence, so a pure
# CDS-density ranking would systematically miss Yan-style regulatory targets.
# The per-gene table reports five signals:
#
#   n_genebody_snps        Focal SNPs in the gene span +/- 2 kb flank.
#                          Primary leaderboard metric. Captures intronic, UTR,
#                          CDS, and immediately-proximal regulatory hits.
#   n_cds_snps             Focal SNPs in the union of CDS exons across all
#                          annotated isoforms. Genes high on both this and
#                          n_genebody_snps are the strongest claim for
#                          "introgressed coding variation."
#   genebody_density       n_genebody_snps / (gene_length_kb + 2*flank_kb).
#                          Normalizes for gene size so a 50 kb gene with 20
#                          SNPs doesn't out-rank a 2 kb gene with 8 SNPs.
#   max_XtXst              Strongest single-SNP signal in the gene+flank
#                          region. Catches genes where one or two very high-
#                          XtXst SNPs would be diluted by density alone.
#   frac_focal_predictor   Fraction of gene-body SNPs whose top predictor is
#                          the locus's nominal env predictor (bio05 for both
#                          chr1 and chr2). Read alongside n_genebody_snps --
#                          a 100% match on 2 SNPs is less informative than
#                          90% on 30.
#
# The CpL3a ranking is exploratory: the inversion, not any gene within it, is
# the unit of inference, and its 11.6 Mb span holds hundreds of genes.
#
# Why +/- 2 kb flank. Standard GFFs annotate protein_coding_gene as TSS-to-TES
# (CDS + UTRs + introns). Core promoters, proximal promoters (~500 bp - 2 kb
# upstream), and immediate 3' regulatory regions are not annotated but are
# biologically attached to the gene. +/- 2 kb captures most of that proximal
# regulatory zone without bleeding into neighboring genes (Culex gene density
# averages ~1 gene per 25 kb, so 2 kb stays inside intergenic space most of
# the time). 5 kb would start double-attributing SNPs between adjacent genes
# and artifactually inflate density. Distal enhancers (10s-100s of kb away)
# are out of reach of this approach -- a known limit of proximity-based
# attribution in non-model insects.
#
# CDS construction. CDS exons are extracted from the GFF, mapped from CDS ->
# mRNA -> gene via the Parent attribute chain, and reduced (union) per gene
# so overlapping isoforms are counted once. cds_length_kb is the total
# length of the gene's CDS bases, not the transcript length.
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")

# ---- Libraries --------------------------------------------------------------

library(tidyverse)
library(rtracklayer)
library(GenomicRanges)


# ---- Configuration ----------------------------------------------------------

BP_FOCAL_FILE <- 'results/baypass/v4/baypass_outliers_clustered.csv'
need_files(BP_FOCAL_FILE)
GFF_FILE      <- deposit_file('data/quinq_on_pip.gff3_polished')

FLANK_BP <- 2000

# Interval settings. filter_to_genes is a gene_id vector to keep, or NULL to
# rank every protein-coding gene in the interval (used for all three).
# focal_predictor is the locus's single nominal env predictor (bio05 for
# CpL1a/CpL2a, pis for CpL3a); focal_category is the broader predictor set
# ('temperature' = bio05/bio07; 'urbanization' = pis/tree_canopy/pd). Both are
# computed per gene because the CpL3a signal is split across tree_canopy and
# pis, so a single-predictor metric undercounts the urbanization signal there.
# The chr2 interval is the 500 kb core around Yan, not the full CpL2a block.
INTERVALS <- list(
  chr1 = list(
    chrom            = CHR1,
    start            = CPL1A_START,
    end              = CPL1A_END,
    filter_to_genes  = NULL,
    focal_predictor  = 'bio05',
    focal_category   = 'temperature',
    out_csv          = 'results/chr1/Table_chr1_cds_snp_attribution.csv'
  ),
  chr2 = list(
    chrom            = CHR2,
    start            = CPL2A_CORE_START,
    end              = CPL2A_CORE_END,
    filter_to_genes  = NULL,
    focal_predictor  = 'bio05',
    focal_category   = 'temperature',
    out_csv          = 'results/chr2/Table_chr2_cds_snp_attribution.csv'
  ),
  chr3 = list(
    chrom            = CHR3,
    start            = CPL3A_START,
    end              = CPL3A_END,
    filter_to_genes  = NULL,
    focal_predictor  = 'pis',
    focal_category   = 'urbanization',
    out_csv          = 'results/chr3/Table_chr3_cds_snp_attribution.csv'
  )
)

# Predictor -> category mapping. Mirrors the predictor_category column produced
# by R/03_baypass/01_baypass_outliers.R. Used to compute
# frac_focal_category per gene in the attribution loop.
PREDICTOR_CATEGORIES <- list(
  temperature  = c('bio05', 'bio07'),
  urbanization = c('pis', 'tree_canopy', 'pd')
)

for (cfg in INTERVALS) dir.create(dirname(cfg$out_csv), recursive = TRUE,
                                  showWarnings = FALSE)


# ---- Load shared inputs -----------------------------------------------------

bp_focal <- read.csv(BP_FOCAL_FILE)
gff_all  <- import.gff3(GFF_FILE)


# =============================================================================
# (1) Build CDS-by-gene index from the GFF
# =============================================================================
# Walk the Parent attribute chain CDS -> mRNA -> gene. Reduce overlapping CDS
# exons within each gene so isoform overlap doesn't double-count. cds_by_gene
# is a named list of GRanges keyed by gene_id.

cds_records  <- gff_all[as.character(gff_all$type) == 'CDS']
mrna_records <- gff_all[as.character(gff_all$type) == 'mRNA']

# CompressedCharacterList parent fields -> character vector of first parent.
first_parent <- function(parent_list) {
  vapply(parent_list, function(x) if (length(x) > 0) as.character(x)[1] else NA_character_,
         character(1))
}

mrna_to_gene <- tibble(
  mrna_id = mcols(mrna_records)$ID,
  gene_id = first_parent(mcols(mrna_records)$Parent)
)

cds_df <- tibble(
  chrom     = as.character(seqnames(cds_records)),
  cds_start = start(cds_records),
  cds_end   = end(cds_records),
  cds_strand= as.character(strand(cds_records)),
  mrna_id   = first_parent(mcols(cds_records)$Parent)
) %>%
  left_join(mrna_to_gene, by = 'mrna_id') %>%
  filter(!is.na(gene_id))

cat(sprintf('GFF: %d CDS records across %d unique genes\n',
            nrow(cds_df), length(unique(cds_df$gene_id))))

# Per-gene reduced CDS GRanges. Returns a list of GRanges keyed by gene_id.
cds_by_gene <- cds_df %>%
  group_split(gene_id) %>%
  setNames(., sapply(., function(d) d$gene_id[1])) %>%
  lapply(function(d) {
    reduce(GRanges(seqnames = d$chrom[1],
                   ranges   = IRanges(start = d$cds_start, end = d$cds_end)))
  })


# =============================================================================
# (2) Family classification
# =============================================================================
# Same regex set as 01_cpl1a_gene_density.R, with transporter added. Used in console output to surface obvious
# candidate gene families without scrolling through descriptions; included as
# a column in the saved CSV.

flag_family <- function(desc) {
  desc_l <- tolower(ifelse(is.na(desc), '', desc))
  case_when(
    str_detect(desc_l, 'cytochrome p450|p450|\\bcyp\\b')                 ~ 'P450',
    str_detect(desc_l, 'glutathione s-transferase|\\bgst\\b')            ~ 'GST',
    str_detect(desc_l, 'carboxyl?esterase|\\bcoe\\b|\\beste(r)?ase\\b')  ~ 'esterase',
    str_detect(desc_l, 'abc transporter|atp-binding cassette|multidrug resistance|multispecific organic anion') ~ 'ABC',
    str_detect(desc_l, 'solute carrier|\\bslc\\d|glucose transporter|sugar transporter|sweet\\d') ~ 'transporter',
    str_detect(desc_l, 'heat shock|chaperone|\\bhsp\\b')                 ~ 'HSP',
    str_detect(desc_l, 'sodium channel|voltage-gated|paralytic')         ~ 'Na_channel',
    str_detect(desc_l, 'potassium channel|\\bshaker\\b|\\bslo\\b')       ~ 'K_channel',
    str_detect(desc_l, 'gaba receptor|\\brdl\\b|chloride channel')       ~ 'GABA/Cl',
    str_detect(desc_l, 'acetylcholinesterase|\\bache\\b|\\bace-?\\d?\\b')~ 'AChE',
    str_detect(desc_l, 'odorant binding|\\bobp\\b')                      ~ 'OBP',
    str_detect(desc_l, 'odorant receptor|olfactory receptor|\\bor\\b')   ~ 'OR',
    str_detect(desc_l, 'ionotropic receptor|\\bir\\d')                   ~ 'IR',
    str_detect(desc_l, 'gustatory receptor|\\bgr\\d')                    ~ 'GR',
    str_detect(desc_l, 'cuticle|chitin')                                 ~ 'cuticle',
    str_detect(desc_l, 'vitellogenin|yolk|chorion')                      ~ 'reproduction',
    str_detect(desc_l, 'defensin|cecropin|attacin|antimicrobial')        ~ 'AMP',
    str_detect(desc_l, 'period|timeless|cryptochrome|circadian|clock')   ~ 'circadian',
    TRUE                                                                 ~ 'other'
  )
}


# =============================================================================
# (3) Per-locus attribution function
# =============================================================================
# Builds the per-gene table for one locus. The body is the same regardless
# of locus; the only differences are the optional gene-id filter and the env
# predictor used for frac_focal_predictor.

attribute_locus <- function(name, chrom, start, end,
                            filter_to_genes, focal_predictor, focal_category,
                            out_csv) {

  cat(sprintf('\n=== %s: %s:%.2f-%.2f Mb ===\n',
              name, chrom, start / 1e6, end / 1e6))

  # ---- Subset focal SNPs to the locus interval ----
  locus_snps <- bp_focal %>%
    filter(chromosome == chrom,
           position   >= start,
           position   <= end)
  cat(sprintf('  Focal SNPs in interval: %d\n', nrow(locus_snps)))
  if (nrow(locus_snps) == 0) {
    warning(sprintf('No focal SNPs in %s locus; skipping.', name))
    return(NULL)
  }

  snps_gr <- GRanges(
    seqnames = locus_snps$chromosome,
    ranges   = IRanges(start = locus_snps$position, end = locus_snps$position),
    XtXst    = locus_snps$XtXst,
    predictor = locus_snps$predictor
  )

  # ---- Subset protein-coding genes to the locus interval ----
  pc_genes <- gff_all[as.character(gff_all$type) == 'protein_coding_gene' &
                      as.character(seqnames(gff_all)) == chrom &
                      start(gff_all) <= end &
                      end(gff_all)   >= start]
  cat(sprintf('  Protein-coding genes in interval: %d\n', length(pc_genes)))

  # Apply gene-id filter if provided.
  if (!is.null(filter_to_genes)) {
    pc_genes <- pc_genes[mcols(pc_genes)$ID %in% filter_to_genes]
    cat(sprintf('  After gene-id filter: %d\n', length(pc_genes)))
    missing <- setdiff(filter_to_genes, mcols(pc_genes)$ID)
    if (length(missing) > 0) {
      warning(sprintf('Genes not found in GFF: %s', paste(missing, collapse = ', ')))
    }
  }

  if (length(pc_genes) == 0) return(NULL)

  # ---- Gene-body + flank GRanges ----
  # Expand each gene's interval by +/- FLANK_BP on each side; clamp at 1.
  genebody_gr <- pc_genes
  start(genebody_gr) <- pmax(1L, start(genebody_gr) - FLANK_BP)
  end(genebody_gr)   <- end(genebody_gr) + FLANK_BP

  # ---- Per-gene loop ----
  per_gene <- lapply(seq_along(pc_genes), function(i) {
    gene_id <- mcols(pc_genes)$ID[i]
    gene_desc <- mcols(pc_genes)$description[i]

    # Gene body + flank overlap
    gb_ovl     <- findOverlaps(snps_gr, genebody_gr[i])
    gb_snps    <- snps_gr[queryHits(gb_ovl)]
    n_gb       <- length(gb_snps)
    gb_size_kb <- width(genebody_gr[i]) / 1e3

    # CDS overlap (only if this gene has CDS records in the index)
    if (gene_id %in% names(cds_by_gene)) {
      gene_cds      <- cds_by_gene[[gene_id]]
      cds_length_kb <- sum(width(gene_cds)) / 1e3
      cds_ovl       <- findOverlaps(snps_gr, gene_cds)
      n_cds         <- length(unique(queryHits(cds_ovl)))
    } else {
      cds_length_kb <- NA_real_
      n_cds         <- NA_integer_
    }

    # Max XtXst within gene+flank, NA if no SNPs
    max_xtx <- if (n_gb > 0) max(mcols(gb_snps)$XtXst, na.rm = TRUE) else NA_real_

    # Fraction of gene-body SNPs on the locus's single nominal predictor.
    frac_focal <- if (n_gb > 0) {
      mean(mcols(gb_snps)$predictor == focal_predictor, na.rm = TRUE)
    } else {
      NA_real_
    }

    # Fraction of gene-body SNPs on any predictor in the locus's focal
    # category (e.g. all temperature predictors, all urbanization predictors).
    # Chr3 in particular has signal split across pis + tree_canopy + pd, so
    # frac_focal_predictor alone (on pis only) undercounts the urbanization
    # signal at the gene level.
    category_predictors <- PREDICTOR_CATEGORIES[[focal_category]]
    frac_focal_cat <- if (n_gb > 0) {
      mean(mcols(gb_snps)$predictor %in% category_predictors, na.rm = TRUE)
    } else {
      NA_real_
    }

    tibble(
      gene_id          = gene_id,
      family           = flag_family(gene_desc),
      gene_start       = start(pc_genes)[i],
      gene_end         = end(pc_genes)[i],
      gene_length_kb   = (end(pc_genes)[i] - start(pc_genes)[i] + 1) / 1e3,
      cds_length_kb    = cds_length_kb,
      n_genebody_snps  = n_gb,
      n_cds_snps       = n_cds,
      genebody_density = n_gb / gb_size_kb,
      cds_density      = if (!is.na(cds_length_kb) && cds_length_kb > 0)
                            n_cds / cds_length_kb else NA_real_,
      max_XtXst        = max_xtx,
      frac_focal_predictor = frac_focal,
      frac_focal_category  = frac_focal_cat,
      description      = gene_desc
    )
  }) %>%
    bind_rows() %>%
    arrange(desc(genebody_density), desc(n_genebody_snps))

  # ---- Save ----
  write.csv(per_gene, out_csv, row.names = FALSE)
  cat(sprintf('  Wrote %s (%d genes)\n', out_csv, nrow(per_gene)))

  per_gene
}


# =============================================================================
# (4) Run for chr1, chr2, chr3
# =============================================================================

results <- list()
for (locus_name in names(INTERVALS)) {
  cfg <- INTERVALS[[locus_name]]
  results[[locus_name]] <- attribute_locus(
    name             = locus_name,
    chrom            = cfg$chrom,
    start            = cfg$start,
    end              = cfg$end,
    filter_to_genes  = cfg$filter_to_genes,
    focal_predictor  = cfg$focal_predictor,
    focal_category   = cfg$focal_category,
    out_csv          = cfg$out_csv
  )
}


# =============================================================================
# (5) Console summary
# =============================================================================
# chr1: print every gene in the order returned (sorted by genebody_density),
# then the ranks of the six CYP6 paralogs.
# chr2/chr3: print top 20 by genebody_density, then top 10 by n_cds_snps,
# then any non-'other' family hits anywhere in the leaderboard. Chr3 output
# is exploratory only (the inversion is the unit of inference).

cat('\n\n=== chr1 (all protein-coding genes in CpL1a) ===\n')
print(results$chr1 %>%
        dplyr::select(gene_id, family, gene_length_kb, cds_length_kb,
                      n_genebody_snps, n_cds_snps,
                      genebody_density, cds_density,
                      max_XtXst, frac_focal_predictor, frac_focal_category),
      n = Inf)

# ---- Leaderboard helper (chr2, chr3) ----
# Same three views for chr2 and chr3 to keep output uniform.

print_leaderboard <- function(res, locus_name) {
  cat(sprintf('\n\n=== %s leaderboard ===\n', locus_name))

  cat('\nTop 20 by genebody_density (with >= 1 gene-body focal SNP):\n')
  print(res %>%
          filter(n_genebody_snps >= 1) %>%
          slice_head(n = 20) %>%
          dplyr::select(gene_id, family,
                        n_genebody_snps, n_cds_snps,
                        genebody_density, cds_density,
                        max_XtXst, frac_focal_predictor, frac_focal_category,
                        description),
        n = 20)

  cat(sprintf('\nTop 10 by n_cds_snps (within %s block):\n', locus_name))
  print(res %>%
          filter(!is.na(n_cds_snps), n_cds_snps >= 1) %>%
          arrange(desc(n_cds_snps), desc(genebody_density)) %>%
          slice_head(n = 10) %>%
          dplyr::select(gene_id, family,
                        n_genebody_snps, n_cds_snps,
                        cds_density, max_XtXst,
                        frac_focal_predictor, frac_focal_category,
                        description),
        n = 10)

  cat('\nNon-"other" family hits with >= 1 gene-body focal SNP:\n')
  fam_hits <- res %>%
    filter(n_genebody_snps >= 1, family != 'other') %>%
    arrange(desc(genebody_density))
  if (nrow(fam_hits) == 0) {
    cat('  (none)\n')
  } else {
    print(fam_hits %>%
            dplyr::select(gene_id, family,
                          n_genebody_snps, n_cds_snps,
                          genebody_density, max_XtXst,
                          frac_focal_predictor, frac_focal_category,
                          description),
          n = Inf)
  }
}

print_leaderboard(results$chr1, 'chr1')
chr1_ranked <- results$chr1 %>%
  mutate(rank_density = row_number(),                 # already sorted by genebody_density
         rank_xtx     = rank(-max_XtXst, ties.method = 'min'))
cat(sprintf('\n=== CYP6 paralogs vs the full %d-gene chr1 interval ===\n', nrow(chr1_ranked)))
print(chr1_ranked %>% filter(gene_id %in% CYP6_GENES) %>%
        dplyr::select(gene_id, rank_density, rank_xtx, n_genebody_snps,
                      n_cds_snps, max_XtXst, frac_focal_predictor), n = Inf)
print_leaderboard(results$chr2, 'chr2')
print_leaderboard(results$chr3, 'chr3')

# =============================================================================
# (6) Targeted check: outlier SNPs in/around Vgsc (para), chromosome 2
# =============================================================================

para_gene_gr <- gff_all[as.character(gff_all$type) == 'protein_coding_gene' &
                        mcols(gff_all)$ID == VGSC_GENE]
stopifnot(length(para_gene_gr) == 1)

para_body_gr <- para_gene_gr
start(para_body_gr) <- pmax(1L, start(para_body_gr) - FLANK_BP)
end(para_body_gr)   <- end(para_body_gr) + FLANK_BP

para_cds_gr <- cds_by_gene[[VGSC_GENE]]
if (is.null(para_cds_gr)) para_cds_gr <- GRanges()

all_outliers_gr <- GRanges(
  seqnames  = bp_focal$chromosome,
  ranges    = IRanges(bp_focal$position, width = 1),
  XtXst     = bp_focal$XtXst,
  predictor = bp_focal$predictor
)

hits <- all_outliers_gr[overlapsAny(all_outliers_gr, para_body_gr)]
hits$in_gene <- overlapsAny(hits, para_gene_gr)   # FALSE => in the +/- FLANK_BP flank
hits$in_cds  <- overlapsAny(hits, para_cds_gr)    # TRUE  => coding

cat(sprintf('\n\n=== %s (para / VGSC): %d outlier SNP(s) in gene +/- %d bp ===\n',
            VGSC_GENE, length(hits), FLANK_BP))
if (length(hits) > 0) {
  print(as.data.frame(hits) %>%
          dplyr::select(chrom = seqnames, pos = start,
                        XtXst, predictor, in_gene, in_cds) %>%
          arrange(desc(XtXst)))
}

# =============================================================================
# End of script
# =============================================================================