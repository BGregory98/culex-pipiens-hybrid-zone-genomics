# =============================================================================
# R/03_baypass/02_cluster_inventory.R
# =============================================================================
# Inventories every cluster of BayPass dual-outlier SNPs genome-wide, not only
# the three candidate loci, and measures what each step of the candidate-SNP
# filter removes. Clusters here are chains of dual-outlier (focal) SNPs only,
# unlike 01_baypass_outliers.R, which chains the full coarse-outlier set; a
# threshold sweep over (gap, minimum size) reports how the cluster count
# depends on that choice. Each cluster carries its genes (lifted JHB2020
# annotation, +/- 2 kb) and the coarse-outlier density of its surrounding Mb.
# Figure 3 is drawn by R/07_multi_locus/11_windowed_pool_pca_figure.R.
#
# Manuscript: Figure 3 (data for panel A); Results "A genome scan identifies three candidate
#             regions with distinct architectures" (1,510 dual outliers, 211 or
#             14% removed by the cluster requirement); Methods "Genome--environment
#             association with BayPass" (14% removed, 1.4% of positions retained)
#
# Inputs:  results/baypass/v4/baypass_outliers_clustered.csv  (R/03_baypass/01_baypass_outliers.R)
#          deposit: results/baypass/v4/filtered_baypass_data.csv
#                   (cluster/04_baypass/filter_baypass_results.R)
#          deposit: data/quinq_on_pip.gff3_polished  (Liftoff of JHB2020 onto idCulPipi1)
# Outputs: results/baypass/Table_genome_wide_cluster_inventory.csv
#          results/baypass/Table_focal_filter_behaviour.csv
#          results/baypass/Table_cluster_threshold_sweep.csv
#
# Runtime: about a minute (GFF import).
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")

suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})


# ---- Configuration ----------------------------------------------------------

FOCAL_FILE    <- 'results/baypass/v4/baypass_outliers_clustered.csv'
FILTERED_FILE <- deposit_file('results/baypass/v4/filtered_baypass_data.csv')
GFF_FILE      <- deposit_file('data/quinq_on_pip.gff3_polished')
need_files(FOCAL_FILE)

OUT_INVENTORY <- 'results/baypass/Table_genome_wide_cluster_inventory.csv'
OUT_FILTER    <- 'results/baypass/Table_focal_filter_behaviour.csv'
OUT_SWEEP     <- 'results/baypass/Table_cluster_threshold_sweep.csv'

dir.create('results/baypass', recursive = TRUE, showWarnings = FALSE)
dir.create('figures/final',   recursive = TRUE, showWarnings = FALSE)

# Primary inventory threshold. 50 kb / 4 SNPs is deliberately more permissive
# than the rule used to define the candidate set, so that the table lists
# everything that could be called a cluster. The sweep below reports what
# happens at stricter settings.
PRIMARY_GAP <- 50000L
PRIMARY_MIN <- 4L

SWEEP_GAPS <- c(10000L, 25000L, 50000L, 100000L)
SWEEP_MINS <- c(4L, 5L, 10L, 11L)

# The candidate-set chaining parameters, for the filter-behaviour section.
PAPER_GAP <- 10000L
PAPER_MIN <- 10L          # chains with length > 10, i.e. >= 11, are kept

# Named candidate intervals, generous bounds so shoulder SNPs are attributed to
# their locus rather than appearing as spurious novel clusters. Figure 3 shades
# the locus intervals themselves (named_bands_fig3), not these bounds.
NAMED <- tribble(
  ~locus,  ~chr, ~start,      ~end,
  'CpL1a', 1L,   105.5e6,     114.0e6,
  'CpL2a', 2L,    11.0e6,      16.5e6,
  'CpL3a', 3L,   146.4e6,     164.0e6
)


# ---- Helpers ----------------------------------------------------------------

# Single-linkage chaining on sorted unique positions. Returns a cluster id per
# position, or integer(0) for an empty input.
chain_ids <- function(pos, max_gap) {
  if (length(pos) == 0L) return(integer(0))
  cumsum(c(Inf, diff(pos)) > max_gap)
}

# Enumerate clusters of >= min_size from a data.frame carrying chr/position.
enumerate_clusters <- function(df, max_gap, min_size) {
  out <- lapply(sort(unique(df$chr)), function(k) {
    p <- sort(unique(df$position[df$chr == k]))
    if (length(p) == 0L) return(NULL)
    cid <- chain_ids(p, max_gap)
    keep <- names(which(table(cid) >= min_size))
    if (length(keep) == 0L) return(NULL)
    lapply(keep, function(g) {
      q <- p[as.character(cid) == g]
      data.frame(chr = k, start = min(q), end = max(q), n_snps = length(q))
    }) %>% bind_rows()
  })
  bind_rows(out)
}

assign_locus <- function(chr, start, end) {
  vapply(seq_along(chr), function(i) {
    hit <- NAMED$chr == chr[i] & start[i] >= NAMED$start & end[i] <= NAMED$end
    if (any(hit)) NAMED$locus[which(hit)[1]] else 'unpromoted'
  }, character(1))
}


# =============================================================================
# (1) Load
# =============================================================================

focal <- fread(FOCAL_FILE)
chrom_levels <- sort(unique(focal$chromosome))
focal[, chr := match(chromosome, chrom_levels)]

cat('\n=== Inputs ===\n')
cat(sprintf('Focal (dual-outlier, clustered) SNPs: %d\n', nrow(focal)))

filt <- fread(FILTERED_FILE,
              select = c('chromosome', 'position',
                         'XtX_outlier_type', 'BF_outlier_type'))
filt[, chr := match(chromosome, chrom_levels)]
cat(sprintf('Pre-clustering outlier rows (SNP x predictor): %d\n', nrow(filt)))


# =============================================================================
# (2) What the focal-SNP filter does
# =============================================================================

pos_all <- unique(filt[, .(chr, position)])

# Positions surviving the candidate-set chaining rule, applied as
# 01_baypass_outliers.R applies it: to the full coarse-outlier set.
kept <- rbindlist(lapply(sort(unique(pos_all$chr)), function(k) {
  p <- sort(unique(pos_all$position[pos_all$chr == k]))
  cid <- chain_ids(p, PAPER_GAP)
  good <- names(which(table(cid) > PAPER_MIN))
  if (!length(good)) return(NULL)
  data.table(chr = k, position = p[as.character(cid) %in% good])
}))

dual <- filt[XtX_outlier_type %in% c('minor', 'major') &
             BF_outlier_type  %in% c('minor', 'major')]
dual_pos <- unique(dual[, .(chr, position)])

reproduced <- merge(dual_pos, kept, by = c('chr', 'position'))

filter_tbl <- tibble(
  step = c('Coarse outlier positions entering clustering (XtX OR BF, >= minor)',
           'Positions retained by the >=10-SNP / 10 kb chaining rule',
           'Dual-outlier positions (>= minor for XtX AND BF)',
           'Dual-outlier AND retained by chaining (published focal set)'),
  n = c(nrow(pos_all), nrow(kept), nrow(dual_pos), nrow(reproduced)),
  pct_of_coarse = round(100 * c(nrow(pos_all), nrow(kept),
                                nrow(dual_pos), nrow(reproduced)) /
                          nrow(pos_all), 2))

cat('\n=== What the focal-SNP filter does ===\n')
print(as.data.frame(filter_tbl), row.names = FALSE)

cat(sprintf('\nReproduction check: %d vs published %d -- %s\n',
            nrow(reproduced), nrow(focal),
            ifelse(nrow(reproduced) == nrow(focal), 'EXACT MATCH',
                   'MISMATCH, investigate before using this table')))
cat(sprintf('\nChaining retains %.1f%% of coarse outlier positions and %.1f%% of\n',
            100 * nrow(kept) / nrow(pos_all),
            100 * nrow(reproduced) / nrow(dual_pos)))
cat('dual-outlier positions. The dual-outlier requirement, not the clustering\n')
cat(sprintf('rule, is the selective step (%.1f%% of coarse positions).\n',
            100 * nrow(dual_pos) / nrow(pos_all)))

write.csv(filter_tbl, OUT_FILTER, row.names = FALSE)

# Coarse-outlier density per Mb, for the local-context column.
dens <- pos_all[, .N, by = .(chr, bin = floor(position / 1e6))]
setnames(dens, 'N', 'coarse_outliers_per_Mb')
genome_median_density <- median(dens$coarse_outliers_per_Mb)
cat(sprintf('\nGenome-wide median coarse-outlier density: %d positions/Mb\n',
            genome_median_density))


# =============================================================================
# (3) Threshold sweep -- how many clusters exist, and where
# =============================================================================

sweep <- lapply(SWEEP_GAPS, function(g) {
  lapply(SWEEP_MINS, function(m) {
    cl <- enumerate_clusters(as.data.frame(focal), g, m)
    if (nrow(cl) == 0L) {
      return(tibble(max_gap_kb = g / 1000, min_snps = m, n_clusters = 0L,
                    n_snps_in_clusters = 0L, n_unpromoted_clusters = 0L,
                    n_snps_unpromoted = 0L))
    }
    cl$locus <- assign_locus(cl$chr, cl$start, cl$end)
    tibble(max_gap_kb = g / 1000, min_snps = m,
           n_clusters = nrow(cl),
           n_snps_in_clusters = sum(cl$n_snps),
           n_unpromoted_clusters = sum(cl$locus == 'unpromoted'),
           n_snps_unpromoted = sum(cl$n_snps[cl$locus == 'unpromoted']))
  }) %>% bind_rows()
}) %>% bind_rows()

cat('\n=== Cluster-count sensitivity to threshold (dual-outliers only) ===\n')
print(as.data.frame(sweep), row.names = FALSE)
write.csv(sweep, OUT_SWEEP, row.names = FALSE)

# Per-locus behaviour at the candidate-set chaining setting.
strict <- enumerate_clusters(as.data.frame(focal), PAPER_GAP, PAPER_MIN + 1L)
if (nrow(strict) > 0) {
  strict$locus <- assign_locus(strict$chr, strict$start, strict$end)
  cat(sprintf('\nAt the nominal rule (>= %d dual-outliers within %d kb chains):\n',
              PAPER_MIN + 1L, PAPER_GAP / 1000L))
  print(strict %>% group_by(locus) %>%
          summarize(n_clusters = n(), n_snps = sum(n_snps), .groups = 'drop') %>%
          as.data.frame(), row.names = FALSE)
  absent <- setdiff(NAMED$locus, strict$locus)
  if (length(absent)) {
    cat(sprintf('  NOTE: %s forms no dual-outlier cluster at this threshold.\n',
                paste(absent, collapse = ', ')))
  }
}


# =============================================================================
# (4) Primary inventory
# =============================================================================

inv <- enumerate_clusters(as.data.frame(focal), PRIMARY_GAP, PRIMARY_MIN)
stopifnot(nrow(inv) > 0)

detail <- lapply(seq_len(nrow(inv)), function(i) {
  s <- focal[chr == inv$chr[i] & position >= inv$start[i] & position <= inv$end[i]]
  pt <- sort(table(s$predictor), decreasing = TRUE)
  ct <- sort(table(s$predictor_category), decreasing = TRUE)
  data.frame(
    top_predictor      = names(pt)[1],
    top_predictor_frac = round(as.numeric(pt[1]) / nrow(s), 2),
    top_category       = names(ct)[1],
    top_category_frac  = round(as.numeric(ct[1]) / nrow(s), 2),
    max_XtXst          = round(max(s$XtXst), 1),
    median_XtXst       = round(median(s$XtXst), 1),
    max_BF_dB          = round(max(s$BF.dB.), 1))
}) %>% bind_rows()

inv <- bind_cols(inv, detail) %>%
  mutate(span_kb = round((end - start) / 1e3, 1),
         locus   = assign_locus(chr, start, end)) %>%
  arrange(desc(n_snps))

# Local coarse-outlier density context.
inv$local_coarse_per_Mb <- vapply(seq_len(nrow(inv)), function(i) {
  b <- dens[chr == inv$chr[i] &
            bin >= floor(inv$start[i] / 1e6) &
            bin <= floor(inv$end[i]   / 1e6)]
  if (!nrow(b)) NA_real_ else round(mean(b$coarse_outliers_per_Mb))
}, numeric(1))
inv$density_vs_genome <- round(inv$local_coarse_per_Mb / genome_median_density, 1)

# Gene annotation: genes overlapping each cluster +/- 2 kb.
# Bioconductor packages are attached here, after the dplyr steps above, so
# that their masking of dplyr verbs cannot affect those steps.
suppressMessages({ library(rtracklayer); library(GenomicRanges) })
cat('\nAnnotating clusters with overlapping genes...\n')
gff <- import.gff3(GFF_FILE)
g   <- gff[gff$type == 'protein_coding_gene' | gff$type == 'gene']
gchr <- match(as.character(seqnames(g)), chrom_levels)

gid   <- sub('^gene[:-]', '', as.character(g$ID))
gdesc <- if ('description' %in% colnames(mcols(g))) {
  as.character(g$description)
} else rep(NA_character_, length(g))
# GFF descriptions are URL-escaped; decode the separators that matter.
gdesc <- gsub('%2C', ',', gsub('%20', ' ', gdesc, fixed = TRUE), fixed = TRUE)

annotate_hits <- function(i, field, n_max) {
  hit <- which(gchr == inv$chr[i] &
               end(g) >= inv$start[i] - 2000 &
               start(g) <= inv$end[i] + 2000)
  if (!length(hit)) return('')
  v <- unique(na.omit(field[hit]))
  v <- v[nzchar(v)]
  if (!length(v)) return('')
  paste(head(v, n_max), collapse = '; ')
}

inv$n_genes   <- vapply(seq_len(nrow(inv)), function(i)
  length(which(gchr == inv$chr[i] &
               end(g) >= inv$start[i] - 2000 &
               start(g) <= inv$end[i] + 2000)), integer(1))
inv$genes     <- vapply(seq_len(nrow(inv)), annotate_hits,
                        character(1), field = gid,   n_max = 12L)
inv$gene_desc <- vapply(seq_len(nrow(inv)), annotate_hits,
                        character(1), field = gdesc, n_max = 8L)

inv <- inv %>%
  mutate(chromosome = chrom_levels[chr],
         cluster_id = sprintf('C%02d', row_number())) %>%
  dplyr::select(cluster_id, locus, chromosome, chr, start, end, span_kb, n_snps,
                top_predictor, top_predictor_frac, top_category, top_category_frac,
                max_XtXst, median_XtXst, max_BF_dB,
                local_coarse_per_Mb, density_vs_genome,
                n_genes, genes, gene_desc)

cat(sprintf('\n=== Primary inventory (>= %d dual-outliers within %d kb chains) ===\n',
            PRIMARY_MIN, PRIMARY_GAP / 1000L))
print(inv %>%
        dplyr::select(cluster_id, locus, chr, start, end, span_kb, n_snps,
                      top_predictor, top_predictor_frac, max_XtXst,
                      density_vs_genome) %>%
        as.data.frame(), row.names = FALSE)

cat('\nSummary by locus assignment:\n')
print(inv %>% group_by(locus) %>%
        summarize(n_clusters = n(), n_snps = sum(n_snps), .groups = 'drop') %>%
        as.data.frame(), row.names = FALSE)

n_in_clusters <- sum(inv$n_snps)
cat(sprintf('\nFocal SNPs in any inventoried cluster: %d of %d (%.1f%%)\n',
            n_in_clusters, nrow(focal), 100 * n_in_clusters / nrow(focal)))
cat(sprintf('Focal SNPs not in any cluster at this threshold: %d (%.1f%%)\n',
            nrow(focal) - n_in_clusters,
            100 * (nrow(focal) - n_in_clusters) / nrow(focal)))

write.csv(inv, OUT_INVENTORY, row.names = FALSE)


# Figure 3 is drawn by R/07_multi_locus/11_windowed_pool_pca_figure.R, which
# sets these candidate SNPs (panel A) over the environment-free windowed PCA
# (panel B).

cat('\n=== Outputs ===\n')
for (f in c(OUT_INVENTORY, OUT_FILTER, OUT_SWEEP)) cat(sprintf('  %s\n', f))
