# =============================================================================
# R/08_loci/03_cyp6_coding_consequences.R
# =============================================================================
# Classifies each BayPass candidate SNP inside the CDS of the six CpL1a CYP6
# paralogs as synonymous, nonsynonymous or nonsense with predictCoding
# (VariantAnnotation) against transcript models built from the lifted
# annotation, and joins the alternate-allele frequency in the three reference
# pools. A CDS SNP can be silent, so this is the test of whether the candidate
# coding variation changes the protein.
#
# Manuscript: Table S6; Results "Cx. quinquefasciatus ancestry and sequencing
#             read depth at a CYP6 cluster rise with summer temperature" (13
#             nonsynonymous); Supplementary Results "Reference structure at the
#             CpL1a candidate SNPs, and why the QUI-versus-PIP contrast is
#             unavailable" (58 CDS SNPs: 13 nonsynonymous, 45 synonymous);
#             Supplementary Methods "Coding consequences, phylogenetic
#             placement, microsynteny and consensus sequences of the CpL1a CYP6
#             paralogs"
#
# Inputs:  results/baypass/v4/baypass_outliers_clustered.csv   (R/03_baypass/01_baypass_outliers.R)
#          deposit: data/quinq_on_pip.gff3_polished
#          deposit: data/reference_genome/cxpip_genome.sanitized.fa (+ .fai)
#          deposit: data/allele_counts_v4.tsv
#          deposit: data/allele_frequencies_v4.csv
# Outputs: results/chr1/Table_chr1_cds_consequences.csv          (per SNP x transcript; Table S6)
#          results/chr1/Table_chr1_cds_consequence_summary.csv   (per paralog counts)
#
# Runtime: ~10 min (TxDb build, allele-count and frequency tables)
# =============================================================================
#
# Design.
#   1. Build a TxDb from the polished lifted annotation. predictCoding() needs
#      transcript structure + CDS phase to place each SNP in frame and on the
#      correct strand. REF/ALT alone is insufficient -- translating a codon
#      needs the two reference bases flanking the SNP -- hence the FASTA.
#   2. Reconstruct the CDS-SNP set of the attribution table by intersecting the
#      BayPass candidate SNPs with the per-gene CDS union of the six paralogs
#      (same construction as 02_candidate_snp_attribution.R). The script prints
#      the per-paralog count so it can be checked against n_cds_snps there.
#   3. predictCoding() returns, per SNP x transcript, the ref/var codon and
#      amino acid and a CONSEQUENCE label.
#   4. Polarization. idCulPipi1 is the Cx. pipiens (background) reference, so REF
#      is the non-introgressed state and ALT (= varAllele) is the candidate
#      QUI-like / introgressed allele. VARAA is therefore the residue encoded by
#      the introgressed allele. FL1/ME1/CAL1 alt-allele frequencies are joined
#      per SNP so the introgressed direction can be confirmed at each site
#      (FL1 high + CAL1 low => ALT is the QUI-like allele; a missense call then
#      reads as "the introgressed allele changes the residue").
#
# Limitation. Pool-seq has no phased haplotypes, so each SNP is classified
# independently. If two candidate SNPs share a codon (codon 503 of the
# CYP6AA7-like paralog), each is scored against the reference codon.
#
# Reference FASTA. predictCoding reads only the positions it needs, so the
# whole-genome FASTA is used. Its chromosome names must match the GFF's
# underscore form (e.g. ENA_OZ004311_OZ004311.1); the seqname guard below stops
# with the FASTA's names if they do not.
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")

library(tidyverse)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)    # cdsBy, extractTranscriptSeqs
library(VariantAnnotation)  # predictCoding
library(Rsamtools)          # FaFile
library(Biostrings)


# ---- Configuration ----------------------------------------------------------

BP_FOCAL   <- 'results/baypass/v4/baypass_outliers_clustered.csv'
need_files(BP_FOCAL)
GFF_FILE   <- deposit_file('data/quinq_on_pip.gff3_polished')
FASTA_FILE <- deposit_file('data/reference_genome/cxpip_genome.sanitized.fa')
COUNTS_TSV <- deposit_file('data/allele_counts_v4.tsv')       # cols: chr,pos,ref,alt,<27 RD:AD>
AF_CSV     <- deposit_file('data/allele_frequencies_v4.csv')  # per-sample ALT-allele freq

dir.create('results/chr1', recursive = TRUE, showWarnings = FALSE)


# ---- (1) TxDb + reference handle --------------------------------------------
# NOTE: on recent Bioconductor, makeTxDbFromGFF moved from GenomicFeatures to
# the txdbmaker package. Use whichever is installed.
make_txdb <- if (requireNamespace('txdbmaker', quietly = TRUE)) {
  txdbmaker::makeTxDbFromGFF
} else {
  GenomicFeatures::makeTxDbFromGFF
}
txdb <- suppressWarnings(make_txdb(GFF_FILE, format = 'gff3'))

# Build the .fai index on first run if absent, then open.
if (!file.exists(paste0(FASTA_FILE, '.fai'))) indexFa(FASTA_FILE)
fa <- FaFile(FASTA_FILE)
open(fa)

# Seqname sanity: TxDb, FASTA, and the SNP tables must share the chr1 label, or
# predictCoding returns zero coding rows silently.
stopifnot(CHR1 %in% seqlevels(txdb))
fa_seqs <- as.character(seqnames(seqinfo(fa)))
if (!CHR1 %in% fa_seqs)
  stop(sprintf(paste0(
    "Chr1 seqname '%s' not found in the FASTA.\n",
    "  FASTA contains names like: %s\n",
    "  -> the sanitized headers differ from the GFF's underscore form; ",
    "rename the FASTA headers to match."),
    CHR1, paste(head(fa_seqs), collapse = ', ')))


# ---- (1b) Frame sanity check ------------------------------------------------
# predictCoding's translation is only correct if the GFF CDS phases are right.
# Translate one CYP6 paralog's CDS from the reference; it should start with M
# and carry no internal stop. If not, the TxDb frame is wrong and every
# consequence call below is untrustworthy -- stop and fix the GFF first.

gff <- import.gff3(GFF_FILE)
first_parent <- function(p) vapply(p, function(x)
  if (length(x)) as.character(x)[1] else NA_character_, character(1))

mrna  <- gff[as.character(gff$type) == 'mRNA']
m2g   <- tibble(mrna_id = mcols(mrna)$ID,
                gene_id = first_parent(mcols(mrna)$Parent))

cds_by_tx <- cdsBy(txdb, by = 'tx', use.names = TRUE)
check_tx  <- m2g %>% filter(gene_id %in% CYP6_GENES) %>% pull(mrna_id) %>%
  intersect(names(cds_by_tx)) %>% head(1)
check_pep <- as.character(translate(
  extractTranscriptSeqs(fa, cds_by_tx[check_tx])[[1]]))
n_internal_stop <- str_count(substr(check_pep, 1, nchar(check_pep) - 1), '\\*')
cat(sprintf('Frame check %s: start=%s  internal_stops=%d  (expect M / 0)\n',
            check_tx, substr(check_pep, 1, 1), n_internal_stop))
if (substr(check_pep, 1, 1) != 'M' || n_internal_stop > 0)
  warning('Frame check failed -- GFF CDS phase is suspect; do not trust calls.')


# ---- (2) Reconstruct the CDS-SNP set at the six paralogs ---------------------
# Per-gene CDS union (reduce across isoforms), restricted to the six paralogs --
# identical construction to 02_candidate_snp_attribution.R.

cdsr   <- gff[as.character(gff$type) == 'CDS']
cds_df <- tibble(chrom   = as.character(seqnames(cdsr)),
                 start   = start(cdsr), end = end(cdsr),
                 mrna_id = first_parent(mcols(cdsr)$Parent)) %>%
  left_join(m2g, by = 'mrna_id') %>%
  filter(gene_id %in% CYP6_GENES)

cds_by_gene <- cds_df %>% group_split(gene_id) %>%
  setNames(sapply(., function(d) d$gene_id[1])) %>%
  lapply(function(d) reduce(GRanges(d$chrom[1], IRanges(d$start, d$end))))

# gene_id -> readable description (CYP6XX-like)
pcg <- gff[as.character(gff$type) == 'protein_coding_gene']
gene_desc <- tibble(gene_id     = mcols(pcg)$ID,
                    description = mcols(pcg)$description) %>%
  filter(gene_id %in% CYP6_GENES)

# Focal SNPs in the CpL1a interval that overlap any paralog CDS.
focal <- read.csv(BP_FOCAL) %>%
  filter(chromosome == CHR1, position >= CPL1A_START, position <= CPL1A_END) %>%
  dplyr::select(chromosome, position, XtXst, predictor)

focal_gr     <- GRanges(focal$chromosome, IRanges(focal$position, width = 1))
paralog_cds  <- unlist(GRangesList(cds_by_gene))
cds_gene_vec <- rep(names(cds_by_gene), lengths(cds_by_gene))  # CDS-exon -> gene

hit <- findOverlaps(focal_gr, paralog_cds)
focal$gene_id <- NA_character_
focal$gene_id[queryHits(hit)] <- cds_gene_vec[subjectHits(hit)]  # paralogs don't overlap

cds_snps <- focal %>% filter(!is.na(gene_id)) %>%
  distinct(chromosome, position, .keep_all = TRUE)

cat(sprintf('\nFocal CDS SNPs at the six paralogs: %d\n', nrow(cds_snps)))
cat('Per-paralog count (should match the attribution table n_cds_snps):\n')
print(cds_snps %>% dplyr::count(gene_id, name = 'n_cds_snps') %>%
        left_join(gene_desc, by = 'gene_id'))


# ---- (3) Attach REF/ALT and donor allele frequencies ------------------------
# REF/ALT from allele_counts_v4.tsv (no header; cols 1-4 = chr,pos,ref,alt).
counts <- fread(COUNTS_TSV, header = FALSE, select = 1:4,
                col.names = c('chromosome', 'position', 'ref', 'alt'))
cds_snps <- cds_snps %>% left_join(counts, by = c('chromosome', 'position'))

# ALT-allele frequency in each reference pool (allele_frequencies_v4.csv stores
# alt-allele frequency per sample).
af <- fread(AF_CSV, select = c('chromosome', 'position', 'FL1', 'ME1', 'CAL1'))
setnames(af, c('FL1', 'ME1', 'CAL1'),
         c('FL1_altfreq', 'ME1_altfreq', 'CAL1_altfreq'))
cds_snps <- cds_snps %>% left_join(af, by = c('chromosome', 'position'))

stopifnot(!any(is.na(cds_snps$alt)))   # every focal SNP must carry an allele


# ---- (4) Predict coding consequences ----------------------------------------
# Query strand '*'; varAllele is the ALT allele on the + strand (predictCoding
# reverse-complements internally for minus-strand transcripts).
query <- GRanges(cds_snps$chromosome,
                 IRanges(cds_snps$position, width = 1),
                 strand = '*')

coding <- predictCoding(query     = query,
                        subject   = txdb,
                        seqSource = fa,
                        varAllele = DNAStringSet(cds_snps$alt))

coding_df <- as.data.frame(coding) %>%
  transmute(
    queryid     = QUERYID,                        # index into `query` / cds_snps
    tx_id       = TXID,
    tx_gene     = GENEID,
    consequence = as.character(CONSEQUENCE),
    refcodon    = as.character(REFCODON),
    varcodon    = as.character(VARCODON),
    refaa       = as.character(REFAA),
    varaa       = as.character(VARAA),
    protein_pos = vapply(PROTEINLOC, function(x) paste(x, collapse = ','),
                         character(1))
  )

# Per SNP x transcript table.
out <- cds_snps %>%
  mutate(queryid = row_number()) %>%
  left_join(coding_df, by = 'queryid') %>%
  left_join(gene_desc, by = 'gene_id')

n_no_call <- out %>% filter(is.na(consequence)) %>%
  distinct(chromosome, position) %>% nrow()
if (n_no_call > 0)
  warning(sprintf('%d CDS SNP(s) returned no coding consequence ',
                  n_no_call), '(annotation/isoform edge -- inspect).')


# ---- (5) Most-severe consequence per SNP, then per-paralog summary ----------
# A SNP can be synonymous in one isoform and missense in another; count it once
# at its most severe call.
sev_map <- c('synonymous' = 1, 'nonsynonymous' = 2,
             'nonsense' = 3, 'frameshift' = 3, 'not translated' = 0)
out$sev <- sev_map[out$consequence]
out$sev[is.na(out$sev)] <- 0

per_snp <- out %>%
  group_by(chromosome, position) %>%
  slice_max(sev, n = 1, with_ties = FALSE) %>%
  ungroup()

per_gene <- per_snp %>%
  group_by(gene_id, description) %>%
  summarise(
    n_cds_snps      = n(),
    n_synonymous    = sum(consequence == 'synonymous',    na.rm = TRUE),
    n_nonsynonymous = sum(consequence == 'nonsynonymous', na.rm = TRUE),
    n_nonsense      = sum(consequence %in% c('nonsense', 'frameshift'),
                          na.rm = TRUE),
    aa_changes      = paste(na.omit(ifelse(
                        consequence == 'nonsynonymous',
                        paste0(refaa, protein_pos, varaa), NA)),
                        collapse = '; '),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_nonsynonymous))


# ---- (6) Write + console summary --------------------------------------------
write.csv(out %>% dplyr::select(
            chromosome, position, gene_id, description, tx_id,
            XtXst, predictor, ref, alt,
            FL1_altfreq, ME1_altfreq, CAL1_altfreq,
            consequence, refcodon, varcodon, refaa, varaa, protein_pos),
          'results/chr1/Table_chr1_cds_consequences.csv', row.names = FALSE)

write.csv(per_gene,
          'results/chr1/Table_chr1_cds_consequence_summary.csv',
          row.names = FALSE)

cat('\n============================================================\n')
cat('CYP6 CDS-SNP coding consequences\n')
cat('============================================================\n')
cat(sprintf('Total focal CDS SNPs (six paralogs): %d\n', nrow(per_snp)))
cat(sprintf('  synonymous:    %d\n', sum(per_snp$consequence == 'synonymous',    na.rm = TRUE)))
cat(sprintf('  nonsynonymous: %d\n', sum(per_snp$consequence == 'nonsynonymous', na.rm = TRUE)))
cat(sprintf('  nonsense:      %d\n', sum(per_snp$consequence %in% c('nonsense','frameshift'), na.rm = TRUE)))
cat(sprintf('Paralogs carrying >=1 nonsynonymous SNP: %d / 6\n',
            sum(per_gene$n_nonsynonymous > 0)))
cat('\nPer-paralog:\n')
print(per_gene)
cat('============================================================\n\n')

# =============================================================================
# End of script
# =============================================================================