# =============================================================================
# R/08_loci/12_cpl3a_chemosensory_enrichment.R
# =============================================================================
# Tests whether the CpL3a inversion is enriched for chemosensory receptor genes
# relative to the rest of the genome, using the curated chemosensory reference
# of 11_chemosensory_gene_reference.R, and compares its yield with the lifted
# JHB2020 annotation. Gene, tandem-array and per-family counts in the inversion
# are compared with the shared 10,000 width-matched random autosomal regions;
# a conditional test randomizes the positions of the inversion's candidate SNPs
# (uniform and circular-shift nulls, 10,000 permutations each).
#
# Manuscript: Supplementary Results "The CpL3a inversion is not enriched for
#             chemosensory receptor genes"; Supplementary Methods "Chemosensory
#             gene content of the CpL3a inversion"; Results "A pre-adapted
#             molestus inversion tracks urbanization independently of
#             temperature" (not enriched)
#
# Inputs:  results/chemosensory/chemosensory_gene_reference.csv   (11_chemosensory_gene_reference.R)
#          results/baypass/v4/baypass_outliers_clustered.csv      (R/03_baypass/01_baypass_outliers.R)
#          results/multilocus/genome_wide_null.rds                (cluster/07_multilocus/01_genome_wide_null.R)
#          results/chr3/Table_chr3_cds_snp_attribution.csv        (02_candidate_snp_attribution.R)
# Outputs: results/chemosensory/Table_CpL3a_chemosensory_genes.csv
#          results/chemosensory/Table_CpL3a_chemosensory_enrichment.csv
#          results/chemosensory/Table_CpL3a_annotation_comparison.csv
#
# Runtime: ~2 min
# =============================================================================
#
# WHY THE CURATED REFERENCE
# -------------------------
# The JHB2020 Liftoff annotation (data/quinq_on_pip.gff3_polished) classifies
# gene families by regex over RefSeq description strings, and RefSeq JHB2020
# under-annotates insect chemosensory receptors -- which is why Menna et al.
# (2026) hand-curated them. The question is whether the block is enriched, or
# whether an 11.6 Mb window simply contains many genes.
#
# THE NULL
# --------
# The SAME 10,000 width-matched random autosomal regions as the other
# locus-level diagnostics (results/multilocus/genome_wide_null.rds, chr3
# configuration: the CpL3a width, rejection-sampled to exclude all three
# candidate intervals). Reusing them rather than resampling makes this test
# directly commensurable with the PCA-FST and per-pool-ancestry permutation
# p-values, with no second sampler.
#
# WHY A CLUSTER-LEVEL STATISTIC TOO
# ---------------------------------
# Chemosensory receptors sit in tandem arrays. A gene-count statistic can be
# driven by one array duplicating locally, which is a single evolutionary
# event, not 20 independent ones. We therefore report gene count AND the
# number of distinct arrays (single-linkage merge at ARRAY_GAP_BP), and test
# both. If only the gene count is significant, the honest reading is "one big
# array", not "enriched for chemoreception".
#
# CAVEAT CARRIED FROM THE REFERENCE TABLE
# ---------------------------------------
# 11_chemosensory_gene_reference.R keeps the LEFTMOST lifted copy per gene
# symbol (recording n_copies). Positional counts here are therefore per-symbol,
# not per-copy, and are conservative wherever a symbol lifted to several
# tandem copies. The script reports how much of the interval this affects.
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")

# ---- Libraries --------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)

set.seed(42)   # conditional permutations in section (4b)


# ---- Configuration ----------------------------------------------------------

CHEM_REF_FILE <- 'results/chemosensory/chemosensory_gene_reference.csv'
BP_FOCAL_FILE <- 'results/baypass/v4/baypass_outliers_clustered.csv'
NULL_RDS_FILE <- 'results/multilocus/genome_wide_null.rds'
JHB_CHR3_FILE <- 'results/chr3/Table_chr3_cds_snp_attribution.csv'
need_files(CHEM_REF_FILE, BP_FOCAL_FILE, NULL_RDS_FILE, JHB_CHR3_FILE)

OUT_GENES   <- 'results/chemosensory/Table_CpL3a_chemosensory_genes.csv'
OUT_ENRICH  <- 'results/chemosensory/Table_CpL3a_chemosensory_enrichment.csv'
OUT_ANNCMP  <- 'results/chemosensory/Table_CpL3a_annotation_comparison.csv'

dir.create('results/chemosensory', recursive = TRUE, showWarnings = FALSE)

# CpL3a is CPL3A_START-CPL3A_END on CHR3 (R/helpers/loci.R). The null regions
# were drawn at this width from the same definition, so the focal region and
# the null regions stay width-matched.

# Gene-body flank for focal-SNP attribution. Matches 02_candidate_snp_attribution.R
# so the two tables are directly comparable.
FLANK_BP <- 2000L

# Single-linkage gap for collapsing tandem arrays into one unit.
# 100 kb is well above Culex intergenic spacing (~25 kb/gene) so genuinely
# adjacent receptors merge, while separate arrays stay separate.
ARRAY_GAP_BP <- 100000L

CHEMO_FAMILIES <- c('IR', 'OR', 'GR', 'OBP')


# ---- Helpers ----------------------------------------------------------------

# One-sided upper-tail empirical p-value with the standard +1 correction, so a
# statistic never returns p = 0 from a finite null.
emp_p_upper <- function(observed, null_vec) {
  null_vec <- null_vec[is.finite(null_vec)]
  (1 + sum(null_vec >= observed)) / (1 + length(null_vec))
}

# Count reference genes overlapping each of many equal-width regions.
# Returns an integer vector, one entry per region.
count_in_regions <- function(gene_gr, regions_df) {
  reg_gr <- GRanges(seqnames = regions_df$chrom,
                    ranges   = IRanges(start = regions_df$start,
                                       end   = regions_df$end))
  countOverlaps(reg_gr, gene_gr, ignore.strand = TRUE)
}

# Collapse a set of gene spans into single-linkage arrays and return the count.
n_arrays <- function(chrom, start, end, gap = ARRAY_GAP_BP) {
  if (length(start) == 0L) return(0L)
  gr <- GRanges(seqnames = chrom, ranges = IRanges(start = start, end = end))
  length(reduce(gr, min.gapwidth = gap, ignore.strand = TRUE))
}

# Array count for each of many regions. Slower than count_in_regions because
# reduce() is per-region, so it is only run on the statistics that need it.
arrays_in_regions <- function(gene_df, regions_df, gap = ARRAY_GAP_BP) {
  gene_gr <- GRanges(seqnames = gene_df$lift_chr,
                     ranges   = IRanges(gene_df$lift_start, gene_df$lift_end))
  reg_gr  <- GRanges(seqnames = regions_df$chrom,
                     ranges   = IRanges(regions_df$start, regions_df$end))
  hits    <- findOverlaps(reg_gr, gene_gr, ignore.strand = TRUE)
  hl      <- split(subjectHits(hits), factor(queryHits(hits),
                                             levels = seq_len(nrow(regions_df))))
  vapply(hl, function(idx) {
    if (length(idx) == 0L) return(0L)
    length(reduce(gene_gr[idx], min.gapwidth = gap, ignore.strand = TRUE))
  }, integer(1), USE.NAMES = FALSE)
}


# =============================================================================
# (1) Load inputs
# =============================================================================

chem_ref <- read.csv(CHEM_REF_FILE, stringsAsFactors = FALSE)
bp_focal <- read.csv(BP_FOCAL_FILE, stringsAsFactors = FALSE)
null_obj <- readRDS(NULL_RDS_FILE)
jhb_chr3 <- read.csv(JHB_CHR3_FILE, stringsAsFactors = FALSE)

cat('\n=== Inputs ===\n')
cat(sprintf('Chemosensory reference genes:      %d (%d lifted onto idCulPipi1)\n',
            nrow(chem_ref), sum(chem_ref$lifted)))
cat(sprintf('v4 focal SNPs genome-wide:         %d\n', nrow(bp_focal)))

# Restrict the reference to lifted, autosomal, family-assigned genes. The
# permutation regions are autosomal, so the observed statistic must be too.
chem <- chem_ref %>%
  filter(lifted,
         gene_family %in% CHEMO_FAMILIES,
         lift_chr %in% c(CHR1, CHR2, CHR3))

cat(sprintf('Autosomal, family-assigned genes:  %d\n', nrow(chem)))
cat('\nGenome-wide family composition of the reference:\n')
print(chem %>% count(gene_family, name = 'n_genes') %>% arrange(desc(n_genes)))

chem_gr <- GRanges(seqnames = chem$lift_chr,
                   ranges   = IRanges(chem$lift_start, chem$lift_end))


# The shared null: 10,000 regions of exactly the inversion's width, already
# rejection-sampled to exclude all three candidate intervals.
null_regions <- null_obj$null$chr3$primary %>%
  as.data.frame() %>%
  dplyr::select(chrom, start, end)

stopifnot(nrow(null_regions) > 0)
# Widths are compared span-inclusive (end - start + 1) on both sides, so the
# check does not trip on a one-base interval-convention difference.
null_w  <- unique(null_regions$end - null_regions$start + 1L)
focal_w <- CPL3A_END - CPL3A_START + 1L
cat(sprintf('\nNull regions: %d, width %s bp (focal width %s bp)\n',
            nrow(null_regions),
            paste(format(null_w, big.mark = ','), collapse = '/'),
            format(focal_w, big.mark = ',')))
if (length(null_w) != 1L || abs(null_w - focal_w) > 1L) {
  warning('Null region width does not match the focal inversion width. ',
          'The permutation p-values below are not width-matched.')
}


# =============================================================================
# (2) Observed chemosensory content of the CpL3a inversion
# =============================================================================

inv_gr <- GRanges(seqnames = CHR3, ranges = IRanges(CPL3A_START, CPL3A_END))

in_inv <- chem %>%
  filter(lift_chr == CHR3,
         lift_end   >= CPL3A_START,
         lift_start <= CPL3A_END) %>%
  arrange(lift_start)

obs_n_genes  <- nrow(in_inv)
obs_n_arrays <- n_arrays(in_inv$lift_chr, in_inv$lift_start, in_inv$lift_end)

cat(sprintf('\n=== Observed: CpL3a inversion (%.2f-%.2f Mb) ===\n',
            CPL3A_START / 1e6, CPL3A_END / 1e6))
cat(sprintf('Chemosensory genes in inversion:   %d\n', obs_n_genes))
cat(sprintf('Distinct tandem arrays (<=%.0f kb): %d\n',
            ARRAY_GAP_BP / 1e3, obs_n_arrays))
cat('\nBy family:\n')
print(in_inv %>% count(gene_family, name = 'n_genes') %>% arrange(desc(n_genes)))

multi_copy <- sum(in_inv$n_copies > 1, na.rm = TRUE)
cat(sprintf('\nSymbols in the interval that lifted to >1 copy: %d',
            multi_copy))
cat(' (counted once each here -- positional counts are conservative)\n')


# =============================================================================
# (3) Attribute v4 focal SNPs to chemosensory gene bodies (+/- 2 kb)
# =============================================================================

focal_chr3 <- bp_focal %>%
  filter(chromosome == CHR3,
         position >= CPL3A_START,
         position <= CPL3A_END)

cat(sprintf('\nv4 focal SNPs inside the inversion: %d\n', nrow(focal_chr3)))

attribute_snps <- function(genes_df, snps_df, flank = FLANK_BP) {
  if (nrow(genes_df) == 0L) return(genes_df)
  g <- GRanges(seqnames = genes_df$lift_chr,
               ranges   = IRanges(pmax(1L, genes_df$lift_start - flank),
                                  genes_df$lift_end + flank))
  s <- GRanges(seqnames = snps_df$chromosome,
               ranges   = IRanges(snps_df$position, snps_df$position))
  hits <- findOverlaps(g, s, ignore.strand = TRUE)

  agg <- tibble(gene_idx = queryHits(hits), snp_idx = subjectHits(hits)) %>%
    mutate(XtXst     = snps_df$XtXst[snp_idx],
           predictor = snps_df$predictor[snp_idx]) %>%
    group_by(gene_idx) %>%
    summarize(n_genebody_snps = n(),
              max_XtXst       = max(XtXst),
              top_predictor   = predictor[which.max(XtXst)],
              .groups = 'drop')

  genes_df %>%
    mutate(gene_idx = row_number()) %>%
    left_join(agg, by = 'gene_idx') %>%
    mutate(n_genebody_snps = replace_na(n_genebody_snps, 0L)) %>%
    dplyr::select(-gene_idx)
}

in_inv <- attribute_snps(in_inv, focal_chr3)
obs_n_with_snp <- sum(in_inv$n_genebody_snps > 0)

cat(sprintf('Chemosensory genes carrying >=1 focal SNP (+/-%d kb): %d of %d\n',
            FLANK_BP / 1000L, obs_n_with_snp, obs_n_genes))


# =============================================================================
# (4) Permutation tests against the shared genome-wide null
# =============================================================================

cat('\n=== Permutation tests (', nrow(null_regions),
    ' width-matched random autosomal regions) ===\n', sep = '')

null_genes  <- count_in_regions(chem_gr, null_regions)
null_arrays <- arrays_in_regions(chem, null_regions)

# Per-family gene counts.
null_by_family <- lapply(CHEMO_FAMILIES, function(fam) {
  fam_gr <- chem_gr[chem$gene_family == fam]
  count_in_regions(fam_gr, null_regions)
})
names(null_by_family) <- CHEMO_FAMILIES

# Genes carrying a focal SNP, tested BETWEEN regions. Reported for
# completeness, but read it with care: focal SNPs are heavily concentrated at
# CpL1a, so most random regions contain zero focal SNPs and this statistic
# mostly measures "does the region have any focal SNPs at all" rather than
# "do its focal SNPs land in chemosensory genes". Section 4b conditions on the
# focal-SNP count and is the test that actually answers the second question.
focal_all_gr <- GRanges(seqnames = bp_focal$chromosome,
                        ranges   = IRanges(bp_focal$position, bp_focal$position))
chem_flank_gr <- GRanges(seqnames = chem$lift_chr,
                         ranges   = IRanges(pmax(1L, chem$lift_start - FLANK_BP),
                                            chem$lift_end + FLANK_BP))
chem_has_focal <- countOverlaps(chem_flank_gr, focal_all_gr,
                                ignore.strand = TRUE) > 0
null_with_snp <- count_in_regions(chem_gr[chem_has_focal], null_regions)

build_row <- function(stat, obs, nullv) {
  tibble(statistic     = stat,
         observed      = obs,
         null_median   = median(nullv),
         null_mean     = round(mean(nullv), 2),
         null_p95      = as.numeric(quantile(nullv, 0.95)),
         null_max      = max(nullv),
         fold_vs_median = ifelse(median(nullv) > 0,
                                 round(obs / median(nullv), 2), NA_real_),
         perm_p        = emp_p_upper(obs, nullv))
}

enrich <- bind_rows(
  build_row('Chemosensory genes',                obs_n_genes,    null_genes),
  build_row('Distinct tandem arrays',            obs_n_arrays,   null_arrays),
  build_row('Genes carrying >=1 focal SNP',      obs_n_with_snp, null_with_snp),
  bind_rows(lapply(CHEMO_FAMILIES, function(fam) {
    build_row(paste0(fam, ' genes'),
              sum(in_inv$gene_family == fam),
              null_by_family[[fam]])
  }))
)

print(as.data.frame(enrich), row.names = FALSE, digits = 3)


# =============================================================================
# (4b) Conditional test: given the inversion's focal SNPs, do they land in
#      chemosensory genes more often than chance?
# =============================================================================
#
# The between-region test in (4) cannot separate "this region has focal SNPs"
# from "this region's focal SNPs hit chemosensory genes". Here we hold the
# observed focal-SNP count inside the inversion fixed and ask only where they
# fall, using two nulls with different assumptions:
#
#   uniform  -- resample all n positions uniformly across the inversion.
#               Simple, but assumes focal SNPs could have landed anywhere,
#               ignoring that they arrive in tight clusters.
#   circular -- shift the observed position vector by a random offset, wrapping
#               at the inversion boundaries. Preserves the observed clustering
#               and spacing, so it is the more conservative of the two.
#
# Neither conditions on where callable/analysable SNPs actually are. Both
# therefore assume the analysable-site density is roughly even across the
# inversion.

N_COND <- 10000L

chem_inv_gr <- GRanges(
  seqnames = in_inv$lift_chr,
  ranges   = IRanges(pmax(1L, in_inv$lift_start - FLANK_BP),
                     in_inv$lift_end + FLANK_BP))

n_focal   <- nrow(focal_chr3)
inv_width <- CPL3A_END - CPL3A_START + 1L

# Number of DISTINCT chemosensory genes hit -- the quantity the manuscript's
# sentence is really about.
genes_hit <- function(pos) {
  gr <- GRanges(seqnames = CHR3, ranges = IRanges(pos, pos))
  sum(countOverlaps(chem_inv_gr, gr, ignore.strand = TRUE) > 0)
}

obs_genes_hit <- genes_hit(focal_chr3$position)

null_unif <- vapply(seq_len(N_COND), function(i)
  genes_hit(sample.int(inv_width, n_focal, replace = TRUE) + CPL3A_START - 1L),
  integer(1))

rel_pos <- focal_chr3$position - CPL3A_START
null_circ <- vapply(seq_len(N_COND), function(i) {
  off <- sample.int(inv_width, 1L)
  genes_hit(((rel_pos + off) %% inv_width) + CPL3A_START)
}, integer(1))

cond <- bind_rows(
  build_row('Chemosensory genes hit (uniform null)',  obs_genes_hit, null_unif),
  build_row('Chemosensory genes hit (circular null)', obs_genes_hit, null_circ)
)

cat(sprintf('\n=== Conditional on the %d focal SNPs inside the inversion ===\n',
            n_focal))
print(as.data.frame(cond), row.names = FALSE, digits = 3)

enrich <- bind_rows(enrich, cond)
write.csv(enrich, OUT_ENRICH, row.names = FALSE)


# =============================================================================
# (5) Annotation comparison: what the manuscript's annotation sees
# =============================================================================
#
# jhb_chr3 is the JHB2020-Liftoff attribution table over CpL3a
# (02_candidate_snp_attribution.R; Table S8), the source of the manuscript's
# "GR2a and Or67d" sentence. Its `family` column is regex over RefSeq
# descriptions. Compare its chemosensory yield to Menna's
# curated annotation over the same interval.

jhb_chemo <- jhb_chr3 %>%
  filter(family %in% CHEMO_FAMILIES)

ann_cmp <- bind_rows(
  tibble(annotation = 'JHB2020 Liftoff (used for Table S8)',
         family     = c(CHEMO_FAMILIES, 'TOTAL'),
         n_genes    = c(vapply(CHEMO_FAMILIES,
                               function(f) sum(jhb_chr3$family == f), integer(1)),
                        nrow(jhb_chemo))),
  tibble(annotation = 'Menna et al. 2026 curated (this script)',
         family     = c(CHEMO_FAMILIES, 'TOTAL'),
         n_genes    = c(vapply(CHEMO_FAMILIES,
                               function(f) sum(in_inv$gene_family == f), integer(1)),
                        obs_n_genes))
) %>%
  pivot_wider(names_from = annotation, values_from = n_genes)

cat('\n=== Annotation comparison over the same interval ===\n')
print(as.data.frame(ann_cmp), row.names = FALSE)
write.csv(ann_cmp, OUT_ANNCMP, row.names = FALSE)

cat('\nThe two genes named in the manuscript, as the JHB2020 annotation sees them:\n')
print(jhb_chemo %>%
        dplyr::select(gene_id, family, gene_start, n_genebody_snps,
                      n_cds_snps, max_XtXst, description) %>%
        as.data.frame(), row.names = FALSE)


# =============================================================================
# (5b) Where the chr3 chemosensory genes actually sit relative to the block
# =============================================================================
#
# A 140-170 Mb window is roughly 2.6x the inversion, so genes in that window
# are not "genes in the inversion"; record the breakdown explicitly.

WIDE_START <- 140e6
WIDE_END   <- 170e6

wide <- chem %>%
  filter(lift_chr == CHR3, lift_start >= WIDE_START, lift_end <= WIDE_END) %>%
  mutate(zone = case_when(
    lift_start >= CPL3A_START & lift_end <= CPL3A_END ~ 'inside inversion',
    lift_start >  CPL3A_END   & lift_start < CPL3A_END + 3e6 ~ 'within 3 Mb distal of breakpoint',
    lift_end   <  CPL3A_START & lift_end   > CPL3A_START - 3e6 ~ 'within 3 Mb proximal of breakpoint',
    TRUE ~ 'elsewhere in 140-170 Mb')) %>%
  arrange(lift_start)

cat('\n=== chr3 chemosensory genes, 140-170 Mb, by position relative to the block ===\n')
print(wide %>% count(zone, name = 'n_genes') %>% arrange(desc(n_genes)) %>%
        as.data.frame(), row.names = FALSE)
cat('\nBy zone and family:\n')
print(wide %>% count(zone, gene_family, name = 'n') %>%
        pivot_wider(names_from = gene_family, values_from = n, values_fill = 0) %>%
        as.data.frame(), row.names = FALSE)

# The distal array is the striking feature; quantify its span and composition.
distal <- wide %>% filter(zone == 'within 3 Mb distal of breakpoint')
if (nrow(distal) > 0) {
  cat(sprintf(paste0('\nDistal array: %d genes spanning %.2f-%.2f Mb ',
                     '(%.2f Mb), starting %.0f kb past the breakpoint.\n'),
              nrow(distal), min(distal$lift_start) / 1e6,
              max(distal$lift_end) / 1e6,
              (max(distal$lift_end) - min(distal$lift_start)) / 1e6,
              (min(distal$lift_start) - CPL3A_END) / 1e3))
  cat('  Composition: ')
  cat(paste(names(table(distal$gene_family)), table(distal$gene_family),
            sep = ' = ', collapse = ', '), '\n')
  cat('  This array is OUTSIDE the block. The manuscript reports LD collapsing\n')
  cat('  to chromosomal baseline immediately past the breakpoints, so it is not\n')
  cat('  carried with the inversion haplotype.\n')
}


# =============================================================================
# (6) Per-gene output table
# =============================================================================

genes_out <- in_inv %>%
  mutate(published_evidence = in_bell | deg_menna | divergent_menna) %>%
  dplyr::select(gene_symbol, gene_family,
                lift_chr, lift_start, lift_end, lift_strand, n_copies,
                in_bell, deg_menna, divergent_menna, menna_annotated,
                published_evidence,
                n_genebody_snps, max_XtXst, top_predictor) %>%
  arrange(lift_start)

write.csv(genes_out, OUT_GENES, row.names = FALSE)

cat(sprintf('\nGenes in the inversion with published functional evidence: %d\n',
            sum(genes_out$published_evidence)))
print(genes_out %>%
        filter(published_evidence) %>%
        dplyr::select(gene_symbol, gene_family, lift_start,
                      in_bell, deg_menna, divergent_menna,
                      n_genebody_snps, max_XtXst, top_predictor) %>%
        as.data.frame(), row.names = FALSE)

cat('\n=== Outputs ===\n')
cat(sprintf('  %s\n', OUT_GENES))
cat(sprintf('  %s\n', OUT_ENRICH))
cat(sprintf('  %s\n', OUT_ANNCMP))
