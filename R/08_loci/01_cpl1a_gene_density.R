# =============================================================================
# R/08_loci/01_cpl1a_gene_density.R
# =============================================================================
# Lists every gene-level feature of the lifted JHB2020 annotation (protein-coding
# genes, pseudogenes, ncRNAs, transposable elements) that overlaps the CpL1a
# interval, flags candidate families by regex on the GFF description, and counts
# the BayPass candidate SNPs inside each feature, with the highest-XtXst SNP and
# its predictor. The table also defines the CpL1a genes that
# 15_kdr_cpl1a_matched_null.R excludes from its gene-level null.
#
# Manuscript: Supplementary Results "A characterized resistance allele outside
#             the candidate set" (through 15_kdr_cpl1a_matched_null.R)
#
# Inputs:  results/baypass/v4/baypass_outliers_clustered.csv   (R/03_baypass/01_baypass_outliers.R)
#          deposit: data/quinq_on_pip.gff3_polished
# Outputs: results/chr1/Table_chr1_gene_intersection.csv
#
# Runtime: ~1 min (GFF import)
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")

# ---- Libraries --------------------------------------------------------------

library(tidyverse)
library(rtracklayer)
library(GenomicRanges)


# ---- Configuration ----------------------------------------------------------

GENE_FLANK <- 0    # SNPs are counted within the feature span itself

BP_FOCAL_FILE <- 'results/baypass/v4/baypass_outliers_clustered.csv'
need_files(BP_FOCAL_FILE)
GFF_FILE      <- deposit_file('data/quinq_on_pip.gff3_polished')

TABLE_OUT <- 'results/chr1/Table_chr1_gene_intersection.csv'
dir.create(dirname(TABLE_OUT), recursive = TRUE, showWarnings = FALSE)


# ---- Load inputs ------------------------------------------------------------

bp_focal <- read.csv(BP_FOCAL_FILE)
gff      <- import.gff3(GFF_FILE)


# ---- Filter to gene-level features (any biotype) ---------------------------
# Every gene-level biotype is kept, not only 'protein_coding_gene', so a
# lncRNA, miRNA, pseudogene or tRNA inside a SNP-dense window is not dropped
# before attribution. Sub-features (mRNA, exon, CDS, UTRs) are excluded because
# they nest inside their parent and would double-count overlaps.

# Feature types in the lifted GFF, genome-wide and in the CpL1a interval.
cat('\n--- Feature types in lifted GFF (genome-wide) ---\n')
print(sort(table(as.character(gff$type)), decreasing = TRUE))

gff_in_interval <- gff[as.character(seqnames(gff)) == CHR1 &
                       start(gff) <= CPL1A_END & end(gff) >= CPL1A_START]
cat(sprintf('\n--- Feature types in chr1 candidate window (%.2f-%.2f Mb) ---\n',
            CPL1A_START / 1e6, CPL1A_END / 1e6))
print(sort(table(as.character(gff_in_interval$type)), decreasing = TRUE))

gene_level_types <- c(
  'protein_coding_gene', 'pseudogene', 'ncRNA_gene', 'gene',
  'lnc_RNA', 'lncRNA', 'ncRNA', 'miRNA', 'tRNA', 'rRNA',
  'snoRNA', 'snRNA', 'antisense_RNA', 'pre_miRNA', 'SRP_RNA',
  'RNase_P_RNA', 'RNase_MRP_RNA', 'telomerase_RNA',
  'transposable_element', 'transposable_element_gene'
)
present_types <- intersect(gene_level_types, unique(as.character(gff$type)))
cat('\n--- Gene-level types retained for annotation ---\n')
cat(paste(' ', present_types, collapse = '\n'), '\n')

genes <- gff[as.character(gff$type) %in% present_types]
cat(sprintf('\nTotal gene-level features genome-wide: %d\n', length(genes)))


# =============================================================================
# (1) Subset focal SNPs and genes to the chr1 candidate interval
# =============================================================================

chr1_snps <- bp_focal %>%
  filter(chromosome == CHR1,
         position   >= CPL1A_START,
         position   <= CPL1A_END) %>%
  mutate(position_Mb = position / 1e6)

# Subset genes whose lifted span overlaps the interval.
genes_in <- genes[seqnames(genes) == CHR1 &
                  start(genes)    <= CPL1A_END   &
                  end(genes)      >= CPL1A_START]

chr1_genes <- tibble(
  gene_id     = mcols(genes_in)$ID,
  gene_type   = as.character(mcols(genes_in)$type),
  gene_desc   = mcols(genes_in)$description,
  gene_start  = start(genes_in),
  gene_end    = end(genes_in),
  gene_strand = as.character(strand(genes_in))
) %>%
  mutate(
    gene_length_kb    = (gene_end - gene_start + 1) / 1e3,
    midpoint_Mb       = (gene_start + gene_end) / 2 / 1e6,
    start_Mb          = gene_start / 1e6,
    end_Mb            = gene_end / 1e6,
    is_protein_coding = gene_type == 'protein_coding_gene'
  )

cat(sprintf('\nchr1 candidate interval: %.2f-%.2f Mb\n',
            CPL1A_START / 1e6, CPL1A_END / 1e6))
cat(sprintf('  Cluster1 focal SNPs in interval:        %d\n', nrow(chr1_snps)))
cat(sprintf('  Gene-level features in interval:        %d\n', nrow(chr1_genes)))
cat(sprintf('    of which protein-coding:              %d\n',
            sum(chr1_genes$is_protein_coding)))
cat(sprintf('    of which non-coding/pseudogene/TE:    %d\n',
            sum(!chr1_genes$is_protein_coding)))
cat('\n--- Biotype breakdown of features in interval ---\n')
print(sort(table(chr1_genes$gene_type), decreasing = TRUE))


# =============================================================================
# (2) Family classification by regex on description
# =============================================================================
#
# These regexes are deliberately broad: they catch any description that
# mentions a candidate-relevant term. A family flag marks a description to
# read, not a functional assignment.

flag_family <- function(desc) {
  desc_l <- tolower(ifelse(is.na(desc), '', desc))
  case_when(
    str_detect(desc_l, 'cytochrome p450|p450|\\bcyp\\b')                 ~ 'P450',
    str_detect(desc_l, 'glutathione s-transferase|\\bgst\\b')            ~ 'GST',
    str_detect(desc_l, 'carboxyl?esterase|\\bcoe\\b|\\beste(r)?ase\\b')  ~ 'esterase',
    str_detect(desc_l, 'abc transporter|atp-binding cassette')           ~ 'ABC',
    str_detect(desc_l, 'heat shock|chaperone|\\bhsp\\b')                 ~ 'HSP',
    str_detect(desc_l, 'sodium channel|voltage-gated|paralytic')         ~ 'Na_channel',
    str_detect(desc_l, 'potassium channel|\\bshaker\\b|\\bslo\\b')       ~ 'K_channel',
    str_detect(desc_l, 'gaba receptor|\\brdl\\b|chloride channel')       ~ 'GABA/Cl',
    str_detect(desc_l, 'acetylcholinesterase|\\bache\\b|\\bace-?\\d?\\b')~ 'AChE',
    str_detect(desc_l, 'odorant binding|\\bobp\\b')                      ~ 'OBP',
    str_detect(desc_l, 'odorant receptor|olfactory receptor|\\bor\\b')   ~ 'OR',
    str_detect(desc_l, 'ionotropic receptor|\\bir\\b')                   ~ 'IR',
    str_detect(desc_l, 'gustatory receptor|\\bgr\\b')                    ~ 'GR',
    str_detect(desc_l, 'cuticle|chitin')                                 ~ 'cuticle',
    str_detect(desc_l, 'vitellogenin|juvenile hormone|ecdysone')         ~ 'reproduction',
    str_detect(desc_l, 'defensin|cecropin|antimicrobial|attacin')        ~ 'AMP',
    str_detect(desc_l, 'period|timeless|\\bcycle\\b|cryptochrome|clock') ~ 'circadian',
    TRUE                                                                  ~ 'other'
  )
}

chr1_genes <- chr1_genes %>%
  mutate(family = flag_family(gene_desc))


# =============================================================================
# (3) Per-gene SNP attribution (within +/- GENE_FLANK)
# =============================================================================

gene_gr <- GRanges(
  seqnames = CHR1,
  ranges   = IRanges(
    start = pmax(1, chr1_genes$gene_start - GENE_FLANK),
    end   =     chr1_genes$gene_end       + GENE_FLANK
  )
)

snp_gr <- GRanges(
  seqnames = chr1_snps$chromosome,
  ranges   = IRanges(start = chr1_snps$position, end = chr1_snps$position)
)

hits <- findOverlaps(gene_gr, snp_gr, ignore.strand = TRUE)

# Per-gene SNP count.
snp_counts <- tibble(gene_idx = queryHits(hits)) %>%
  dplyr::count(gene_idx, name = 'n_focal_snps')

# Per-gene best (highest-XtXst) focal SNP.
best_snp <- tibble(
  gene_idx = queryHits(hits),
  snp_idx  = subjectHits(hits)
) %>%
  mutate(
    snp_position      = chr1_snps$position[snp_idx],
    snp_XtXst         = chr1_snps$XtXst[snp_idx],
    snp_predictor     = chr1_snps$predictor[snp_idx],
    snp_predictor_cat = chr1_snps$predictor_category[snp_idx]
  ) %>%
  group_by(gene_idx) %>%
  slice_max(snp_XtXst, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-snp_idx)

chr1_gene_table <- chr1_genes %>%
  mutate(gene_idx = row_number()) %>%
  left_join(snp_counts, by = 'gene_idx') %>%
  left_join(best_snp,   by = 'gene_idx') %>%
  mutate(n_focal_snps = replace_na(n_focal_snps, 0L)) %>%
  mutate(n_focal_snps_per_kb = n_focal_snps/gene_length_kb) %>%
  dplyr::select(gene_id, gene_type, is_protein_coding, family, gene_desc,
                gene_start, gene_end, gene_strand, gene_length_kb,
                n_focal_snps, n_focal_snps_per_kb, snp_position, snp_XtXst,
                snp_predictor, snp_predictor_cat) %>%
  arrange(desc(n_focal_snps_per_kb), desc(n_focal_snps), desc(snp_XtXst))


# =============================================================================
# (4) Diagnostic summaries
# =============================================================================

cat('\n--- Family breakdown (all features in interval) ---\n')
print(as.data.frame(
  chr1_genes %>% dplyr::count(family, name = 'n_genes') %>% arrange(desc(n_genes))
), row.names = FALSE)

cat('\n--- Predictor-category breakdown across cluster1 focal SNPs ---\n')
print(chr1_snps %>%
        dplyr::count(predictor_category, name = 'n_snps') %>%
        mutate(percent = round(100 * n_snps / sum(n_snps), 1)))

cat('\n--- Predictor-level breakdown across cluster1 focal SNPs ---\n')
print(chr1_snps %>%
        dplyr::count(predictor, name = 'n_snps') %>%
        arrange(desc(n_snps)) %>%
        mutate(percent = round(100 * n_snps / sum(n_snps), 1)))

cat('\n--- Top features by n_focal_snps (within +/- ',
    GENE_FLANK / 1e3, ' kb) ---\n', sep = '')
print(chr1_gene_table %>%
        filter(n_focal_snps > 0) %>%
        head(20),
      n = 20)

cat('\n--- Non-"other" family genes in interval ---\n')
print(chr1_gene_table %>%
        filter(family != 'other') %>%
        dplyr::select(gene_id, gene_type, family, gene_desc,
                      gene_start, gene_end,
                      n_focal_snps, snp_XtXst, snp_predictor),
      n = 100)

# ---- Non-protein-coding features --------------------------------------------
# Surface every non-coding/pseudogene/TE feature in the interval, regardless
# of family classification (most non-coding features will fall into 'other'
# under the description-regex flag_family() rules).

cat('\n--- All non-protein-coding features in interval ---\n')
nc_features <- chr1_gene_table %>%
  dplyr::filter(!is_protein_coding) %>%
  dplyr::select(gene_id, gene_type, family, gene_desc,
                gene_start, gene_end, gene_length_kb,
                n_focal_snps, snp_XtXst, snp_predictor)
if (nrow(nc_features) == 0) {
  cat('  (no non-coding, pseudogene, or TE features in interval)\n')
} else {
  print(nc_features, n = Inf)
}

cat('\n--- Non-protein-coding features that overlap focal SNPs ---\n')
nc_with_snps <- nc_features %>% dplyr::filter(n_focal_snps > 0)
if (nrow(nc_with_snps) == 0) {
  cat('  (none -- all focal-SNP overlaps in interval are with protein-coding features)\n')
} else {
  print(nc_with_snps, n = Inf)
}


# =============================================================================
# (5) Write output
# =============================================================================

write.csv(chr1_gene_table, TABLE_OUT, row.names = FALSE)
cat(sprintf('\nWrote %s\n', TABLE_OUT))
