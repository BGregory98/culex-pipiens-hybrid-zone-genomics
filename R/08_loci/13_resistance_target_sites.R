# =============================================================================
# R/08_loci/13_resistance_target_sites.R
# =============================================================================
# Reads the coding variation at canonical insecticide-resistance target and
# detoxification genes directly from the v4 allele counts, independently of
# the BayPass genome scan. Every coding SNP in the target CDS is classified
# with predictCoding on the longest-CDS (canonical) transcript; nonsynonymous
# variants get per-pool frequencies (alt/depth, missing below 10 reads) and OLS
# fits on bio05 and percent impervious surface. The Vgsc candidate is
# reconciled to housefly (Musca domestica Vssc1) numbering by global alignment
# to confirm kdr L1014F.
#
# Manuscript: Table S3; Results "A genome scan identifies three candidate
#             regions with distinct architectures" (63-gene panel, kdr);
#             Supplementary Results "A characterized resistance allele outside
#             the candidate set"; Supplementary Methods "Targeted analysis of
#             insecticide-resistance target sites"
#
# Inputs:  data/site_coords_with_spatial_data.csv          (R/00_setup/01_site_predictors.R)
#          results/resistance/Q94615_Musca_Vssc1.fasta     (UniProt Q94615, Musca domestica Vssc1)
#          deposit: data/quinq_on_pip.gff3_polished
#          deposit: data/reference_genome/cxpip_genome.sanitized.fa (+ .fai)
#          deposit: data/allele_counts_v4.tsv            no header; chr, pos, ref, alt,
#                                                        then 27 pools as RD:AD
#          deposit: data/allele_frequencies_v4.csv       header supplies the pool column order
# Outputs: results/resistance/Table_resistance_target_coding_variants.csv   (per variant)
#          results/resistance/Table_resistance_gene_summary.csv            (Table S3)
#          results/resistance/vgsc_coding_variant_counts.tsv               (read by 14_kdr_frequency.R)
#          results/resistance/Table_vgsc_kdr_confirmation.csv
#          results/resistance/txdb_quinq_on_pip.sqlite   (TxDb cache, rebuilt if absent)
#
# Runtime: ~15 min (TxDb build on first run; one pass over the allele counts)
# =============================================================================
#
# WHY THIS IS A DIFFERENT QUESTION FROM THE GENOME SCAN
# -----------------------------------------------------
# The genome scan asks "is this locus environmentally structured?", and a
# population can be fixed for kdr and show no environmental association at all.
# A fixed resistance allele is a strong operational result and an invisible
# BayPass result. So this script does not filter on outlier status. It reads
# every coding variant in the target genes and reports its frequency.
#
# WHY WE DO NOT HARDCODE THE CANONICAL CODONS
# -------------------------------------------
# The canonical substitutions are numbered in other species' conventions --
# Vgsc L1014F in housefly numbering, Ace-1 G119S in Torpedo numbering, Rdl
# A296S in Drosophila numbering. Transferring those numbers onto a Liftoff-
# projected Culex annotation requires an alignment to the reference protein,
# and getting it silently wrong produces a confidently wrong answer at exactly
# the site people care about. Instead we enumerate ALL nonsynonymous variants
# in the target genes with their protein positions and per-pool frequencies,
# then flag those whose residue change matches a canonical substitution
# (L->F, G->S, A->S, R->L, R->Q) as candidates for confirmation. Position-level
# confirmation against the reference protein is a separate, deliberate step
# (section 5, Vgsc only). The enumeration is a superset of the canonical-site
# question and is independent of numbering conventions.
#
# TARGETS
#   Vgsc / para          pyrethroid target (kdr)
#   Ace-1                organophosphate/carbamate target (G119S)
#   Rdl / GABA-gated Cl  cyclodiene/fipronil target (A296S)
#   Ester locus          carboxylesterases, the classic Culex OP-resistance
#                        amplification (esterase B1 / FE4)
#   GSTs                 metabolic detoxification
#   CpL1a CYP6 paralogs  included so CYP6BZ2 R213L/R213Q (Aardema et al. 2022,
#                        ~75% combined frequency across the global complex) is
#                        tested inside CpL1a
#   Na channel 60E       NOT a resistance target; a labelled negative control
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")

suppressMessages({
  library(data.table)
  library(dplyr)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(VariantAnnotation)
  library(Rsamtools)
  library(Biostrings)
})
if (!requireNamespace('pwalign', quietly = TRUE))
  stop("Package 'pwalign' (Bioconductor) is required for the kdr alignment in section 5.")


# ---- Configuration ----------------------------------------------------------

ENV_CSV     <- 'data/site_coords_with_spatial_data.csv'
MUSCA_FASTA <- 'results/resistance/Q94615_Musca_Vssc1.fasta'
need_files(ENV_CSV, MUSCA_FASTA)

GFF_FILE   <- deposit_file('data/quinq_on_pip.gff3_polished')
FASTA_FILE <- deposit_file('data/reference_genome/cxpip_genome.sanitized.fa')
COUNTS_TSV <- deposit_file('data/allele_counts_v4.tsv')
AF_HEADER  <- deposit_file('data/allele_frequencies_v4.csv')
TXDB_CACHE <- 'results/resistance/txdb_quinq_on_pip.sqlite'

OUT_VARIANTS <- 'results/resistance/Table_resistance_target_coding_variants.csv'
OUT_SUMMARY  <- 'results/resistance/Table_resistance_gene_summary.csv'
OUT_VGSC_CNT <- 'results/resistance/vgsc_coding_variant_counts.tsv'
OUT_KDR      <- 'results/resistance/Table_vgsc_kdr_confirmation.csv'

dir.create('results/resistance', recursive = TRUE, showWarnings = FALSE)

CHUNK_ROWS <- 2e6L        # rows per pass over the 11.6M-row counts file
MIN_DEPTH  <- 10L         # per-pool depth below which a frequency is set NA

# Named gene sets. Regex is matched against the GFF description field.
TARGET_PATTERNS <- list(
  'Vgsc/para (kdr)'        = 'sodium channel protein para',
  'Ace-1'                  = 'acetylcholinesterase',
  'Rdl / GABA-gated Cl'    = 'gamma-aminobutyric acid receptor subunit',
  'Ester locus'            = 'esterase B|esterase FE4|carboxylesterase',
  'GST'                    = 'glutathione S-transferase',
  'Na channel 60E (control)' = 'sodium channel protein 60E'
)

# Residue changes that match a canonical resistance substitution. Direction is
# not enforced (the reference may carry either state), so both are listed.
CANONICAL <- c('L>F','F>L','G>S','S>G','A>S','S>A','R>L','L>R','R>Q','Q>R')


# ---- (1) Annotation, targets, reference -------------------------------------

cat('\n=== Loading annotation ===\n')
gff  <- import.gff3(GFF_FILE)
gns  <- gff[gff$type %in% c('gene','protein_coding_gene')]
gid  <- sub('^gene[:-]', '', as.character(gns$ID))
gdsc <- gsub('%2C', ',', gsub('%20', ' ', as.character(gns$description), fixed = TRUE),
             fixed = TRUE)

targets <- rbindlist(lapply(names(TARGET_PATTERNS), function(nm) {
  i <- grep(TARGET_PATTERNS[[nm]], gdsc, ignore.case = TRUE)
  if (!length(i)) return(NULL)
  data.table(target_class = nm, gene_id = gid[i],
             chr = as.character(seqnames(gns))[i],
             gstart = start(gns)[i], gend = end(gns)[i], desc = gdsc[i])
}))
targets <- rbind(targets,
  data.table(target_class = 'CpL1a CYP6 paralog',
             gene_id = CYP6_GENES,
             chr = as.character(seqnames(gns))[match(CYP6_GENES, gid)],
             gstart = start(gns)[match(CYP6_GENES, gid)],
             gend   = end(gns)[match(CYP6_GENES, gid)],
             desc   = gdsc[match(CYP6_GENES, gid)]))
targets <- unique(targets, by = 'gene_id')

cat(sprintf('Target genes: %d\n', nrow(targets)))
print(targets[, .N, by = target_class])

# TxDb: expensive to build from a 55 MB GFF, so cache it.
make_txdb <- if (requireNamespace('txdbmaker', quietly = TRUE))
  txdbmaker::makeTxDbFromGFF else GenomicFeatures::makeTxDbFromGFF
if (file.exists(TXDB_CACHE)) {
  cat('Loading cached TxDb\n')
  txdb <- loadDb(TXDB_CACHE)
} else {
  cat('Building TxDb from GFF (slow, one time)...\n')
  txdb <- suppressWarnings(make_txdb(GFF_FILE, format = 'gff3'))
  saveDb(txdb, TXDB_CACHE)
}

fa <- FaFile(FASTA_FILE)
if (!file.exists(paste0(FASTA_FILE, '.fai'))) indexFa(fa)

cds_by_gene <- cdsBy(txdb, by = 'gene')
names(cds_by_gene) <- sub('^gene[:-]', '', names(cds_by_gene))
have <- intersect(targets$gene_id, names(cds_by_gene))
cat(sprintf('Target genes with CDS in TxDb: %d of %d\n', length(have), nrow(targets)))

cds_gr <- unlist(reduce(cds_by_gene[have]))
cds_gr$gene_id <- rep(have, lengths(reduce(cds_by_gene[have])))
cat(sprintf('Total target CDS bases: %s\n',
            format(sum(width(cds_gr)), big.mark = ',')))


# ---- (2) Pull allele counts at target CDS positions -------------------------

pool_names <- strsplit(readLines(AF_HEADER, n = 1), ',')[[1]][-(1:3)]
n_pools    <- length(pool_names)
cat(sprintf('\nPools: %d (ME1/FL1/CAL1 at %s)\n', n_pools,
            paste(match(c('ME1','FL1','CAL1'), pool_names), collapse = '/')))

cat('Scanning allele-count file in chunks...\n')
keep <- list(); skip <- 0L; total <- 0L
repeat {
  ch <- tryCatch(fread(COUNTS_TSV, skip = skip, nrows = CHUNK_ROWS,
                       header = FALSE, showProgress = FALSE),
                 error = function(e) NULL)
  if (is.null(ch) || nrow(ch) == 0L) break
  total <- total + nrow(ch)
  setnames(ch, 1:4, c('chr','pos','ref','alt'))
  hit <- overlapsAny(GRanges(ch$chr, IRanges(ch$pos, ch$pos)), cds_gr)
  if (any(hit)) keep[[length(keep) + 1L]] <- ch[hit]
  skip <- skip + CHUNK_ROWS
  cat(sprintf('  %s rows scanned, %s retained\r', format(total, big.mark = ','),
              format(sum(vapply(keep, nrow, integer(1))), big.mark = ',')))
  if (nrow(ch) < CHUNK_ROWS) break
}
cat('\n')
dat <- rbindlist(keep)
setnames(dat, 5:(4 + n_pools), pool_names)
cat(sprintf('SNPs in target CDS: %s (from %s total)\n',
            format(nrow(dat), big.mark = ','), format(total, big.mark = ',')))
stopifnot(nrow(dat) > 0)


# ---- (3) Predict coding consequences ----------------------------------------

cat('\n=== Predicting coding consequences ===\n')
query <- GRanges(dat$chr, IRanges(dat$pos, dat$pos),
                 strand = '*',
                 REF = DNAStringSet(dat$ref),
                 varAllele = DNAStringSet(dat$alt))
names(query) <- paste0(dat$chr, ':', dat$pos)

coding <- predictCoding(query = query, subject = txdb, seqSource = fa,
                        varAllele = query$varAllele)
cat(sprintf('Coding rows returned: %d\n', length(coding)))

cdt <- as.data.table(mcols(coding))
cdt[, `:=`(chr = as.character(seqnames(coding)), pos = start(coding))]
cdt[, gene_id := sub('^gene[:-]', '', as.character(GENEID))]
cdt <- cdt[gene_id %in% have]

# CANONICAL TRANSCRIPT. These genes carry many Liftoff isoforms of differing
# length -- para alone has 51, spanning 2,039-2,148 aa -- and a residue index
# means nothing without saying which isoform it indexes. Taking the min or
# first PROTEINLOC across isoforms yields a number that is not comparable
# between genes and cannot be checked against any published numbering, which is
# exactly how a confidently wrong site call gets made. We therefore designate
# the longest-CDS transcript per gene as canonical and report positions in it.
cds_by_tx <- cdsBy(txdb, by = 'tx', use.names = TRUE)
tx_len    <- sum(width(cds_by_tx))
tx_meta   <- as.data.table(transcripts(txdb, columns = c('tx_id','tx_name','gene_id')))
tx_meta[, gene_id := sub('^gene[:-]', '', as.character(gene_id))]
tx_meta[, cds_len := tx_len[match(tx_name, names(tx_len))]]
canon <- tx_meta[!is.na(cds_len) & gene_id %in% have][order(gene_id, -cds_len),
                 .(canon_tx = tx_name[1], canon_aa = floor(cds_len[1] / 3) - 1),
                 by = gene_id]
cat(sprintf('Canonical transcripts designated for %d genes\n', nrow(canon)))

cdt[, tx_name := as.character(tx_meta$tx_name[match(TXID, tx_meta$tx_id)])]
cdt <- merge(cdt, canon, by = 'gene_id', all.x = TRUE)
cdt[, aa_change := paste0(as.character(REFAA), '>', as.character(VARAA))]

# Keep the canonical-transcript row per SNP; fall back to the first isoform
# only when the SNP is not coding in the canonical transcript, and flag it.
cdt[, is_canon := tx_name == canon_tx]
snp_cons <- cdt[order(-is_canon)][
  , .(consequence  = CONSEQUENCE[1],
      protein_pos  = unlist(PROTEINLOC)[1],
      aa_change    = aa_change[1],
      canon_tx     = canon_tx[1],
      canon_aa_len = canon_aa[1],
      on_canonical = is_canon[1],
      n_isoforms   = .N,
      n_iso_calls  = uniqueN(CONSEQUENCE)),
  by = .(chr, pos, gene_id)]

cat(sprintf('SNPs placed on the canonical transcript: %d of %d\n',
            sum(snp_cons$on_canonical), nrow(snp_cons)))

# Peptide context around each nonsynonymous site, so a call can be checked
# against a published motif without trusting our residue index.
cat('Extracting peptide context...\n')
canon_seqs <- extractTranscriptSeqs(fa, cds_by_tx[unique(snp_cons$canon_tx)])
canon_prot <- suppressWarnings(translate(canon_seqs, if.fuzzy.codon = 'solve'))
snp_cons[, peptide_context := {
  p <- canon_prot[[canon_tx]]
  if (is.na(protein_pos) || protein_pos > length(p)) NA_character_
  else as.character(subseq(p, max(1, protein_pos - 8),
                              min(length(p), protein_pos + 8)))
}, by = seq_len(nrow(snp_cons))]

cat('\nConsequence breakdown across target CDS SNPs:\n')
print(snp_cons[, .N, by = consequence][order(-N)])


# ---- (4) Per-pool frequencies and environmental association -----------------

parse_ad <- function(x) {
  s  <- tstrsplit(x, ':', fixed = TRUE)
  rd <- as.numeric(s[[1]]); ad <- as.numeric(s[[2]])
  dp <- rd + ad
  list(freq = fifelse(dp >= MIN_DEPTH, ad / dp, NA_real_), dp = dp)
}

freq <- as.data.table(lapply(dat[, ..pool_names], function(c) parse_ad(c)$freq))
setnames(freq, pool_names)
freq[, `:=`(chr = dat$chr, pos = dat$pos)]

env    <- fread(ENV_CSV)
fields <- intersect(pool_names, env$site)
cat(sprintf('\nField pools with environmental data: %d\n', length(fields)))
env    <- env[match(fields, site)]

ns <- merge(snp_cons[consequence == 'nonsynonymous'], freq, by = c('chr','pos'))
cat(sprintf('Nonsynonymous variants in target genes: %d\n', nrow(ns)))

fit_env <- function(y, x) {
  ok <- is.finite(y) & is.finite(x)
  if (sum(ok) < 8 || sd(y[ok]) == 0) return(c(NA_real_, NA_real_))
  m <- summary(lm(y[ok] ~ x[ok]))$coefficients
  if (nrow(m) < 2) return(c(NA_real_, NA_real_))
  c(m[2,1], m[2,4])
}

fm <- as.matrix(ns[, ..fields])
stats <- rbindlist(lapply(seq_len(nrow(ns)), function(i) {
  y  <- as.numeric(fm[i, ])
  b5 <- fit_env(y, env$bio05); pi <- fit_env(y, env$pis)
  data.table(mean_field_freq = round(mean(y, na.rm = TRUE), 3),
             min_field_freq  = round(min(y,  na.rm = TRUE), 3),
             max_field_freq  = round(max(y,  na.rm = TRUE), 3),
             bio05_slope = round(b5[1], 4), bio05_p = signif(b5[2], 3),
             pis_slope   = round(pi[1], 5), pis_p   = signif(pi[2], 3))
}))

out <- cbind(
  ns[, .(target_class = targets$target_class[match(gene_id, targets$gene_id)],
         gene_id, desc = targets$desc[match(gene_id, targets$gene_id)],
         chr, pos, canon_tx, canon_aa_len, protein_pos, on_canonical,
         aa_change, peptide_context, n_isoforms, n_iso_calls)],
  stats,
  ns[, .(ME1 = round(ME1, 3), FL1 = round(FL1, 3), CAL1 = round(CAL1, 3))])
out[, canonical_residue_match := aa_change %in% CANONICAL]
out <- out[order(target_class, gene_id, protein_pos)]

write.csv(out, OUT_VARIANTS, row.names = FALSE)

cat('\n=== Nonsynonymous variants per target class ===\n')
print(out[, .(n_variants = .N,
              n_genes = uniqueN(gene_id),
              n_canonical_residue_match = sum(canonical_residue_match),
              n_env_assoc_p05 = sum(bio05_p < 0.05 | pis_p < 0.05, na.rm = TRUE)),
          by = target_class][order(-n_variants)])

cat('\n=== Variants whose residue change matches a canonical substitution ===\n')
hits <- out[canonical_residue_match == TRUE][order(-mean_field_freq)]
if (nrow(hits)) {
  print(head(hits[, .(target_class, gene_id, pos, protein_pos, aa_change,
                      mean_field_freq, min_field_freq, max_field_freq,
                      bio05_p, pis_p, FL1, ME1, CAL1)], 25), row.names = FALSE)
} else cat('  none\n')

cat('\n=== Highest-frequency nonsynonymous variants overall ===\n')
print(head(out[order(-mean_field_freq),
               .(target_class, gene_id, protein_pos, aa_change,
                 mean_field_freq, bio05_p, pis_p, FL1, ME1, CAL1)], 20),
      row.names = FALSE)

gene_sum <- out[, .(n_nonsyn = .N,
                    max_mean_freq = max(mean_field_freq, na.rm = TRUE),
                    n_canonical   = sum(canonical_residue_match),
                    min_bio05_p   = suppressWarnings(min(bio05_p, na.rm = TRUE)),
                    min_pis_p     = suppressWarnings(min(pis_p,   na.rm = TRUE))),
                by = .(target_class, gene_id, desc)][order(target_class, -n_nonsyn)]
write.csv(gene_sum, OUT_SUMMARY, row.names = FALSE)

# Raw RD:AD counts at the Vgsc nonsynonymous variants, in the allele-count
# table's layout (chr, pos, ref, alt, 27 pools); 14_kdr_frequency.R reads the
# kdr site from here rather than from the full counts file.
vgsc_sites  <- out[gene_id == VGSC_GENE, .(chr, pos)]
vgsc_counts <- dat[vgsc_sites, on = .(chr, pos)][order(pos)]
fwrite(vgsc_counts, OUT_VGSC_CNT, sep = '\t')


# ---- (5) Confirm the kdr site against the housefly reference ----------------
#
# Turns the residue-match flag at Vgsc into an actual site call. kdr is numbered
# L1014F in Musca domestica Vssc1 (UniProt Q94615); our annotation's canonical
# para transcript uses its own numbering, so the two are reconciled by global
# alignment.

KDR_TX <- canon$canon_tx[canon$gene_id == VGSC_GENE]
stopifnot(length(KDR_TX) == 1L)

cat('\n=== Confirming the Vgsc kdr site against Musca domestica Q94615 ===\n')
mus <- readAAStringSet(MUSCA_FASTA)[[1]]
cxp <- suppressWarnings(translate(
  extractTranscriptSeqs(fa, cds_by_tx[KDR_TX]), if.fuzzy.codon = 'solve'))[[1]]

al <- pwalign::pairwiseAlignment(cxp, mus, substitutionMatrix = 'BLOSUM62',
                                 gapOpening = 10, gapExtension = 0.5,
                                 type = 'global')
cv <- strsplit(as.character(pwalign::alignedPattern(al)), '')[[1]]
sv <- strsplit(as.character(pwalign::alignedSubject(al)), '')[[1]]
ci <- cumsum(cv != '-'); si <- cumsum(sv != '-')

kdr_rows <- out[gene_id == VGSC_GENE]
kdr_rows[, musca_pos := vapply(protein_pos, function(pp) {
  k <- which(ci == pp & cv != '-')[1]
  if (is.na(k)) NA_integer_ else as.integer(si[k])
}, integer(1))]
kdr_rows[, musca_call := ifelse(is.na(musca_pos), NA_character_,
                                paste0(substr(aa_change, 1, 1), musca_pos,
                                       substr(aa_change, 3, 3)))]

cat(sprintf('Global alignment: %.1f%% identity, score %d\n',
            pwalign::pid(al), round(pwalign::score(al))))
print(kdr_rows[, .(protein_pos, aa_change, musca_pos, musca_call,
                   peptide_context, mean_field_freq, FL1, ME1, CAL1)],
      row.names = FALSE)

kdr <- kdr_rows[musca_call == 'L1014F']
if (nrow(kdr) == 1L) {
  cat(sprintf(paste0('\nCONFIRMED kdr L1014F at %s:%s\n',
                     '  field mean %.3f (range %.3f-%.3f) | ',
                     'FL1 %.3f  ME1 %.3f  CAL1 %.3f\n',
                     '  bio05 P = %.3g | impervious surface P = %.3g\n'),
              kdr$chr, format(kdr$pos, big.mark = ','),
              kdr$mean_field_freq, kdr$min_field_freq, kdr$max_field_freq,
              kdr$FL1, kdr$ME1, kdr$CAL1, kdr$bio05_p, kdr$pis_p))
  cat('  Not a BayPass focal SNP: it segregates widely but without an\n')
  cat('  environmental cline, so a GEA scan cannot recover it.\n')
} else {
  cat('\nNo variant mapped to Musca L1014 - do not report a kdr call.\n')
}
write.csv(kdr_rows, OUT_KDR, row.names = FALSE)

cat('\n=== Outputs ===\n')
for (f in c(OUT_VARIANTS, OUT_SUMMARY, OUT_VGSC_CNT, OUT_KDR)) cat(sprintf('  %s\n', f))
cat('\nNOTE: canonical_residue_match flags a matching AMINO-ACID CHANGE, not a\n')
cat('confirmed canonical site. Only Vgsc is reconciled to published numbering\n')
cat('by alignment in step 5. Ace-1 G119S and Rdl A296S residue matches remain\n')
cat('candidates, not calls, until the same step is done for them.\n')
