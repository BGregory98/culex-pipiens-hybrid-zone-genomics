# =============================================================================
# R/08_loci/11_chemosensory_gene_reference.R
# =============================================================================
# Builds one reference table of chemosensory genes from Bell et al. (2024,
# Heredity) and Menna et al. (2026, Genomics), with each gene's coordinates on
# the idCulPipi1 (Cx. pipiens) assembly taken from Menna's manually curated
# chemosensory annotation lifted onto idCulPipi1.
#
# Manuscript: Supplementary Methods "Chemosensory gene content of the CpL3a
#             inversion" (486 gene symbols; the reference used by
#             12_cpl3a_chemosensory_enrichment.R)
#
# Inputs:  data/menna_on_pip.gff3_polished, the curated chemosensory annotation
#          of Menna et al. (2026) on idCulPipi1. It belongs to that publication
#          and is not in this study's data deposit; to rerun, obtain it from
#          Menna et al. (2026) and place it at that path in the deposit
#          (deposit_file). The output table is in the repository, so
#          12_cpl3a_chemosensory_enrichment.R runs without it, and run_all.R
#          skips this script when the file is absent.
# Outputs: results/chemosensory/chemosensory_gene_reference.csv
#
# Runtime: ~1 min (GFF import)
# =============================================================================
#
# Four sources are merged:
#   (i)   Bell et al. Table 3 -- 36 ORs/OBPs in host-preference-associated
#         RDA windows. Coordinates are from JHB2020 and are hardcoded below
#         from the published table for traceability.
#   (ii)  Menna et al. Tables 4-7 -- genes that are differentially expressed
#         in the antennae of avian- vs human-seeking Cx. pipiens. Identified
#         by gene name; coordinates pulled from the lifted GFF3.
#   (iii) Menna et al. Supplementary Data S7 -- genes that are sequence-
#         divergent (phylogenetically) between bioforms in the antennal
#         transcriptome.
#   (iv)  Menna et al. full curated chemosensory annotation -- every gene in
#         the lifted GFF3 whose `locus` attribute starts with Or/Ir/Obp/Gr
#         (plus Orco). This includes Menna-annotated genes that are NOT in
#         her published DEG or divergence tables, which is needed to represent
#         the chemosensory arrays on chromosome 3 near CpL3a: several IRs
#         there (Ir174-Ir186, Ir75l, etc.) are annotated but functionally
#         uncharacterized.
#
# Source flags written into the output:
#   in_bell           -- TRUE if in Bell Table 3
#   deg_menna         -- TRUE if in Menna Tables 4-7 (published DEG)
#   divergent_menna   -- TRUE if in Menna Supplementary Data S7 (sequence
#                        divergence)
#   menna_annotated   -- TRUE if present in Menna's lifted curated annotation
#                        (true for nearly every row -- distinguishes "exists
#                        in annotation" from "Bell-only with no Menna entry")
#
# Plus `lifted` (TRUE if lifted-assembly coordinates were recovered) and
# `gene_family` (OBP/OR/IR/GR/other).
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")

# ---- Libraries --------------------------------------------------------------

library(tidyverse)
library(rtracklayer)
library(GenomicRanges)


# ---- Configuration ----------------------------------------------------------

GFF_FILE <- deposit_file('data/menna_on_pip.gff3_polished')
OUT_FILE <- 'results/chemosensory/chemosensory_gene_reference.csv'
dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

# JHB2020 chromosome NCBI accessions, retained for traceability of Bell
# Table 3 source-table coordinates only. Not used by the matcher.
JHB_CHR_MAP <- c(
  NC_051861.1 = 'chr1',
  NC_051862.1 = 'chr2',
  NC_051863.1 = 'chr3'
)

# Symbol prefixes that mark chemosensory gene families. Apollo placeholder
# names (locus11164, locus15857) are rejected by the `locus[0-9]+` regex,
# since those have no curated annotation behind them.
CHEMO_PREFIX_RE  <- '^(Orco|Or[0-9]|Obp[0-9]|Ir[0-9]|Gr[0-9])'
PLACEHOLDER_RE   <- '^locus[0-9]+$'


# =============================================================================
# (1) Bell et al. (2024) Table 3 -- ORs and OBPs in host-preference RDA windows
# =============================================================================

bell_table3 <- tribble(
  ~gene_symbol, ~jhb_chr,        ~jhb_start, ~jhb_end,    ~strand,
  'Obp92',      'NC_051861.1',   3852588,    3853514,     '+',
  'Obp2',       'NC_051861.1',   51920776,   51921939,    '-',
  'Or202',      'NC_051862.1',   209436027,  209437270,   '-',
  'Or140',      'NC_051862.1',   209437564,  209438850,   '-',
  'Or201',      'NC_051862.1',   209440193,  209441503,   '-',
  'Or142',      'NC_051862.1',   209441711,  209442967,   '-',
  'Or143',      'NC_051862.1',   209444125,  209445381,   '-',
  'Or144',      'NC_051862.1',   209445675,  209446878,   '+',
  'Obp82',      'NC_051862.1',   49806497,   49807372,    '-',
  'Obp81',      'NC_051862.1',   49809142,   49810005,    '+',
  'Obp78',      'NC_051862.1',   49810242,   49811120,    '+',
  'Obp77',      'NC_051862.1',   49811322,   49812188,    '+',
  'Obp76',      'NC_051862.1',   49813384,   49814238,    '+',
  'Obp80',      'NC_051862.1',   49814428,   49815162,    '+',
  'Obp79',      'NC_051862.1',   49819306,   49820466,    '+',
  'Obp83',      'NC_051862.1',   49820632,   49821319,    '+',
  'Obp120',     'NC_051862.1',   49821446,   49822252,    '+',
  'Or22',       'NC_051862.1',   51249865,   51251328,    '+',
  'Or64',       'NC_051862.1',   58393834,   58395729,    '-',
  'Or65',       'NC_051862.1',   58402332,   58408334,    '-',
  'Obp22',      'NC_051862.1',   79148170,   79148879,    '+',
  'Obp21',      'NC_051862.1',   79151997,   79152727,    '-',
  'Or57',       'NC_051862.1',   95233480,   95235109,    '+',
  'Or56',       'NC_051862.1',   95244528,   95250839,    '+',
  'Or55',       'NC_051862.1',   95255042,   95256585,    '-',
  'Or54',       'NC_051862.1',   95260174,   95261760,    '-',
  'Or53',       'NC_051862.1',   95261771,   95264830,    '-',
  'Or52',       'NC_051862.1',   95267126,   95268875,    '-',
  'Or51',       'NC_051862.1',   95292384,   95293992,    '+',
  'Or176',      'NC_051862.1',   95294779,   95296245,    '+',
  'Obp11',      'NC_051863.1',   143710283,  143715125,   '+',
  'Or208',      'NC_051863.1',   143859543,  143860894,   '-',
  'Obp14',      'NC_051863.1',   46700831,   46713009,    '-',
  'Or40',       'NC_051863.1',   46861523,   46863363,    '-',
  'Orco',       'NC_051863.1',   46915625,   46949854,    '-',
  'Obp8',       'NC_051863.1',   47120103,   47120842,    '-'
) %>%
  mutate(
    in_bell  = TRUE,
    jhb_chr_label = JHB_CHR_MAP[jhb_chr]
  )


# =============================================================================
# (2) Menna et al. (2026) Tables 4-7 -- antennally differentially expressed
# =============================================================================

menna_deg_irs <- c(
  'Ir40a','Ir41a','Ir41i','Ir41j','Ir41k','Ir41l','Ir41o','Ir41p',
  'Ir75c','Ir75i','Ir75m','Ir75n','Ir75o','Ir75q','Ir75r',
  'Ir76b','Ir87a','Ir105','Ir171','Ir181','Ir183','Ir31a1'
)

menna_deg_grs <- c('Gr19','Gr23','Gr37','Gr39a','Gr68')

menna_deg_ors <- c(
  'Or1','Or2','Or3','Or4','Or6','Or22','Or25','Or27','Or32','Or35','Or37',
  'Or40','Or41','Or46','Or50','Or53','Or63','Or65','Or66','Or76','Or78',
  'Or80','Or93','Or97','Or105a','Or114','Or115','Or129','Or136','Or160',
  'Or164','Or165','Or167','Or168','Or171','Or173','Or174','Or177','Or181',
  'Or192','Or195','Or196','Or197','Or205'
)

menna_deg_obps <- c(
  'Obp1','Obp2','Obp4','Obp9','Obp10','Obp12','Obp15','Obp17','Obp18a',
  'Obp18b','Obp19','Obp20','Obp24','Obp29','Obp51','Obp63','Obp65','Obp104',
  'Obp118','Obp125'
)

menna_deg_all <- c(menna_deg_irs, menna_deg_grs, menna_deg_ors, menna_deg_obps)

menna_deg <- tibble(
  gene_symbol = menna_deg_all,
  deg_menna   = TRUE
)


# =============================================================================
# (3) Menna et al. (2026) Supplementary Data S7 -- "fully divergent" OBPs/ORs
# =============================================================================

menna_div_obps <- c('Obp17','Obp19','Obp20','Obp24','Obp104')
menna_div_ors  <- c('Or11','Or36','Or119','Or209','Orco','Or194','Or93')

menna_divergent <- tibble(
  gene_symbol     = c(menna_div_obps, menna_div_ors),
  divergent_menna = TRUE
)


# =============================================================================
# (4) Lifted GFF3 -- recover idCulPipi1 coordinates and full chemosensory
#     annotation from Menna's curated set
# =============================================================================

gff   <- import.gff3(GFF_FILE)
genes <- gff[gff$type == 'gene']

if (!'locus' %in% colnames(mcols(genes))) {
  stop("Lifted GFF gene records do not carry a 'locus' attribute. ",
       "Verify that ", GFF_FILE, " was lifted from Menna's curated GFF and ",
       "not from the JHB2020 RefSeq GFF.")
}
n_locus_populated <- sum(!is.na(genes$locus) & nzchar(genes$locus))
cat(sprintf('Loaded %d gene records from %s (%d with locus= attribute).\n',
            length(genes), GFF_FILE, n_locus_populated))

# One-row-per-symbol coordinate lookup. With Liftoff -copies, a symbol may
# appear on multiple lifted records (tandem paralogs). We keep the leftmost
# occurrence and record n_copies for transparency.
locus_lookup <- tibble(
  gene_symbol = as.character(genes$locus),
  lift_chr    = as.character(seqnames(genes)),
  lift_start  = start(genes),
  lift_end    = end(genes),
  lift_strand = as.character(strand(genes))
) %>%
  filter(!is.na(gene_symbol), nzchar(gene_symbol)) %>%
  group_by(gene_symbol) %>%
  mutate(n_copies = n()) %>%
  slice_min(lift_start, n = 1, with_ties = FALSE) %>%
  ungroup()

# Chemosensory-prefixed symbols from Menna's curated annotation. This is the
# source of the `menna_annotated` flag and the source of new rows beyond
# Bell+Menna's published tables.
menna_annotated_symbols <- locus_lookup %>%
  filter(grepl(CHEMO_PREFIX_RE, gene_symbol),
         !grepl(PLACEHOLDER_RE, gene_symbol)) %>%
  pull(gene_symbol) %>%
  unique()

cat(sprintf('Menna-annotated chemosensory genes (Or/Ir/Obp/Gr/Orco): %d\n',
            length(menna_annotated_symbols)))

menna_annotated_tbl <- tibble(
  gene_symbol     = menna_annotated_symbols,
  menna_annotated = TRUE
)


# =============================================================================
# (5) Merge all four sources into a single reference table
# =============================================================================

all_symbols <- unique(c(bell_table3$gene_symbol,
                        menna_deg$gene_symbol,
                        menna_divergent$gene_symbol,
                        menna_annotated_symbols))

ref <- tibble(gene_symbol = all_symbols) %>%
  left_join(locus_lookup,        by = 'gene_symbol') %>%
  left_join(bell_table3 %>% dplyr::select(gene_symbol, in_bell,
                                          jhb_chr_label, jhb_start, jhb_end),
            by = 'gene_symbol') %>%
  left_join(menna_deg,           by = 'gene_symbol') %>%
  left_join(menna_divergent,     by = 'gene_symbol') %>%
  left_join(menna_annotated_tbl, by = 'gene_symbol') %>%
  mutate(
    across(c(in_bell, deg_menna, divergent_menna, menna_annotated),
           ~ tidyr::replace_na(., FALSE)),
    lifted = !is.na(lift_chr),
    gene_family = case_when(
      grepl('^Orco$', gene_symbol, ignore.case = TRUE) ~ 'OR',
      grepl('^Obp',   gene_symbol)                     ~ 'OBP',
      grepl('^Or',    gene_symbol)                     ~ 'OR',
      grepl('^Ir',    gene_symbol)                     ~ 'IR',
      grepl('^Gr',    gene_symbol)                     ~ 'GR',
      TRUE                                             ~ 'other'
    )
  ) %>%
  arrange(gene_family, gene_symbol)


# =============================================================================
# (6) Diagnostic summary
# =============================================================================

cat('\n=== Chemosensory gene reference ===\n')
cat(sprintf('Total unique gene symbols:                       %d\n', nrow(ref)))
cat(sprintf('  in Bell Table 3:                               %d\n',
            sum(ref$in_bell)))
cat(sprintf('  Menna DEG (any family):                        %d\n',
            sum(ref$deg_menna)))
cat(sprintf('  Menna fully divergent (sequence):              %d\n',
            sum(ref$divergent_menna)))
cat(sprintf('  Menna curated chemosensory annotation:         %d\n',
            sum(ref$menna_annotated)))
cat(sprintf('Successfully lifted onto idCulPipi1:             %d (%.1f%%)\n',
            sum(ref$lifted),
            100 * mean(ref$lifted)))

unlifted_published <- ref %>%
  filter(!lifted, in_bell | deg_menna | divergent_menna) %>%
  pull(gene_symbol)

if (length(unlifted_published) > 0) {
  cat('\nUnlifted gene symbols from published tables (',
      length(unlifted_published), ' total):\n', sep = '')
  cat(paste(' ', unlifted_published, collapse = '\n'), '\n')
} else {
  cat('\nAll Bell+Menna published-table gene symbols lifted successfully.\n')
}

cat('\nLift breakdown by gene family:\n')
print(ref %>%
        group_by(gene_family) %>%
        summarize(total  = n(),
                  lifted = sum(lifted),
                  in_bell_or_menna_published = sum(in_bell | deg_menna | divergent_menna),
                  .groups = 'drop'))

# Chromosome 3 genes around CpL3a (140-170 Mb).
cat('\nGenes on chr3 140-170 Mb (around CpL3a):\n')
ref %>%
  filter(lifted,
         lift_chr == CHR3,
         lift_start > 140e6,
         lift_end   < 170e6) %>%
  arrange(lift_start) %>%
  dplyr::select(gene_symbol, gene_family, lift_start, lift_end,
                in_bell, deg_menna, divergent_menna, menna_annotated) %>%
  print(n = Inf)


# =============================================================================
# (7) Write output
# =============================================================================

write.csv(ref, OUT_FILE, row.names = FALSE)
cat(sprintf('\nWrote %s\n', OUT_FILE))