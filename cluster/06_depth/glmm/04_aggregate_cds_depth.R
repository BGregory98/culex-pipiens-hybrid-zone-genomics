#!/usr/bin/env Rscript
# cluster/06_depth/glmm/04_aggregate_cds_depth.R
# ==============================================
#
# Step 4 of the callable-depth pipeline (cluster/README.md step 8; Methods
# "Callable-position read-depth analysis"). Sums callable depth within each
# gene's CDS per pool and the genome-wide callable depth per pool. Submitted
# through 04_aggregate_cds_depth.slurm (field pools) and
# 04b_aggregate_cds_depth_all27.slurm (all 27 pools).
#
# Inputs:
#   --callable-cds      callable_cds_depth.tsv.gz from step 3
#                       (chrom, pos, gene_id, then per-pool depths)
#   --callable-all      callable_depth_v4.tsv.gz from step 2
#                       (chrom, pos, then per-pool depths) - needed for the
#                       genome-wide-callable-depth offset denominator
#   --sample-categories sample_categories.csv
#                       used to classify pools as field or reference
#   --min-positions     minimum callable-position count per gene to retain
#                       (default 50)
#   --pool-set          'field' (default) or 'all'. See below.
#   --out-dir           where to write outputs
#
# Outputs (in --out-dir):
#   cds_depth_per_gene_pool.tsv      gene_id, pool, total_depth, n_callable_pos
#                                    Long format. One row per (gene, pool).
#                                    Filtered to genes with >= min_positions
#                                    callable positions and to field pools only.
#   genome_callable_depth_per_pool.tsv   pool, genome_callable_depth
#                                    One row per pool. Sum of depths across all
#                                    callable positions genome-wide. This is the
#                                    offset denominator in the GLMM.
#
# These two files are the inputs of R/06_depth/01_callable_depth_glmm.R
# (repo data/glmm_v4/).
#
# --pool-set
# ----------
# The GLMM analyses the 24 field pools only, so the default 'field' writes
# the GLMM inputs above.
#
# 'all' additionally emits the three reference pools (QUI/FL1, PIP/ME1,
# MOL/CAL1), writing to *_all27.tsv so the GLMM inputs above are never
# overwritten. This exists for the reference-pool depth control at the CpL1a
# CYP6 cluster: whether the cluster reads as depth-elevated in a pool that
# matches the assembly (PIP) tells you whether elevation is a property of the
# sample or of the reference, which the field pools alone cannot separate.
#
# Step 2 already writes all 27 pool columns, so this script alone, run on
# the existing step-3 output, produces the control.
#
# The callable set is defined by min-cov 10 across ALL 27 pools including the
# references, so field and reference depths are measured on exactly the same
# positions. It also means positions where QUI drops out are already excluded
# genome-wide and at the cluster alike, so this control speaks to copy number
# and not to dropout.
#
# Performance
# -----------
# data.table::fread is the I/O bottleneck. For ~5 GB callable_cds_depth.tsv.gz
# (rough order of magnitude), expect ~5-15 min wall time. The aggregation
# itself is fast.
#
# Usage
# -----
#     Rscript 04_aggregate_cds_depth.R \
#         --callable-cds callable_cds_depth.tsv.gz \
#         --callable-all callable_depth_v4.tsv.gz \
#         --sample-categories sample_categories.csv \
#         --min-positions 50 \
#         --out-dir glmm_aggregation_v4/
#
# and for the reference-pool control:
#
#     Rscript 04_aggregate_cds_depth.R \
#         ... same arguments ... \
#         --pool-set all

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# ---- Argument parsing -------------------------------------------------------

option_list <- list(
  make_option('--callable-cds', type = 'character', default = NULL,
              help = 'callable_cds_depth.tsv.gz from step 3', metavar = 'PATH'),
  make_option('--callable-all', type = 'character', default = NULL,
              help = 'callable_depth_v4.tsv.gz from step 2', metavar = 'PATH'),
  make_option('--sample-categories', type = 'character', default = NULL,
              help = 'sample_categories.csv', metavar = 'PATH'),
  make_option('--min-positions', type = 'integer', default = 50,
              help = 'minimum callable-position count per gene to retain [default %default]'),
  make_option('--pool-set', type = 'character', default = 'field',
              help = "'field' for the 24 field pools (GLMM inputs, default) or 'all' to include the three reference pools [default %default]"),
  make_option('--out-dir', type = 'character', default = NULL,
              help = 'output directory', metavar = 'PATH')
)
opt <- parse_args(OptionParser(option_list = option_list))

for (req in c('callable-cds', 'callable-all', 'sample-categories', 'out-dir')) {
  if (is.null(opt[[req]])) {
    stop(sprintf('--%s is required', req))
  }
}

if (!opt[['pool-set']] %in% c('field', 'all')) {
  stop(sprintf("--pool-set must be 'field' or 'all', got '%s'", opt[['pool-set']]))
}
pool_set <- opt[['pool-set']]

dir.create(opt[['out-dir']], showWarnings = FALSE, recursive = TRUE)


# ---- Sample categories: classify pools --------------------------------------

cat(sprintf('[INFO] loading sample categories from %s\n', opt[['sample-categories']]))
sample_info <- fread(opt[['sample-categories']])
stopifnot('sample' %in% names(sample_info), 'sample_type' %in% names(sample_info))

field_pools <- sort(sample_info$sample[sample_info$sample_type == 'field'])
cat(sprintf('[INFO] %d field pools: %s\n', length(field_pools),
            paste(field_pools, collapse = ', ')))

# Reference pools are everything not typed 'field'. Note that CAL1 (MOL) is
# typed 'colony_p', NOT 'reference_*', so matching on 'reference' would silently
# drop the MOL reference -- the one whose panel carries the CpL3a result. Define
# the set by exclusion, and intersect with the columns actually present, because
# sample_categories.csv lists replicate pools (FL2/FL3, ME2/ME3, Evanston,
# Northfield) that are not in the 27-pool v4 mpileup.
nonfield_pools <- sort(setdiff(sample_info$sample, field_pools))
cat(sprintf('[INFO] %d non-field pools in metadata: %s\n', length(nonfield_pools),
            paste(nonfield_pools, collapse = ', ')))
cat(sprintf("[INFO] --pool-set %s\n", pool_set))


# ---- Load callable CDS depth, restrict to field pools, aggregate -----------

cat(sprintf('[INFO] reading %s\n', opt[['callable-cds']]))
cds_dt <- fread(opt[['callable-cds']])

# Schema check: chrom, pos, gene_id, then >= 27 pool columns.
stopifnot(names(cds_dt)[1:3] == c('chrom', 'pos', 'gene_id'))
all_pools <- names(cds_dt)[4:ncol(cds_dt)]
missing_field <- setdiff(field_pools, all_pools)
if (length(missing_field) > 0) {
  stop(sprintf('field pools missing from callable CDS depth file: %s',
               paste(missing_field, collapse = ', ')))
}

cat(sprintf('[INFO] callable CDS rows: %s\n', format(nrow(cds_dt), big.mark = ',')))
cat(sprintf('[INFO] unique gene_ids:   %s\n', format(uniqueN(cds_dt$gene_id), big.mark = ',')))

# Which pools to carry through. 'field' reproduces the original outputs exactly.
ref_present <- intersect(nonfield_pools, all_pools)
if (pool_set == 'all') {
  if (length(ref_present) == 0) {
    stop('--pool-set all requested but no non-field pool columns are present in ',
         'the callable CDS file. Check that step 2 ran with the 27-pool locations ',
         'file rather than a field-only one.')
  }
  pools_use <- c(field_pools, ref_present)
  cat(sprintf('[INFO] carrying %d pools: %d field + %d reference (%s)\n',
              length(pools_use), length(field_pools), length(ref_present),
              paste(ref_present, collapse = ', ')))
  unused <- setdiff(nonfield_pools, ref_present)
  if (length(unused) > 0) {
    cat(sprintf('[INFO] non-field pools in metadata but absent from this mpileup: %s\n',
                paste(unused, collapse = ', ')))
  }
} else {
  pools_use <- field_pools
  cat(sprintf('[INFO] carrying %d field pools; %d reference column(s) present but excluded (%s)\n',
              length(pools_use), length(ref_present),
              if (length(ref_present)) paste(ref_present, collapse = ', ') else 'none'))
}

# Reshape to long: gene_id, pool, depth.
cds_long <- melt(
  cds_dt,
  id.vars       = c('chrom', 'pos', 'gene_id'),
  measure.vars  = pools_use,
  variable.name = 'pool',
  value.name    = 'depth',
  variable.factor = FALSE
)

# Aggregate to gene x pool.
cat('[INFO] aggregating to gene x pool depth sums\n')
cds_agg <- cds_long[, .(
    total_depth    = sum(depth),
    n_callable_pos = .N
  ), by = .(gene_id, pool)]

# Apply minimum-positions filter at the gene level (some gene has fewer than
# min_positions callable positions => drop from analysis entirely, not just for
# certain pools).
cat(sprintf('[INFO] filtering to genes with >= %d callable positions\n',
            opt[['min-positions']]))
pos_per_gene <- cds_agg[, .(n_callable_pos = first(n_callable_pos)), by = gene_id]
genes_keep   <- pos_per_gene$gene_id[pos_per_gene$n_callable_pos >= opt[['min-positions']]]

cat(sprintf('[INFO] %s genes pass min-positions filter (of %s total)\n',
            format(length(genes_keep), big.mark = ','),
            format(nrow(pos_per_gene), big.mark = ',')))

cds_agg_filt <- cds_agg[gene_id %in% genes_keep]


# ---- Load callable-all, sum genome-wide per pool ---------------------------

cat(sprintf('[INFO] reading %s for genome-wide callable depth offset\n',
            opt[['callable-all']]))
all_dt <- fread(opt[['callable-all']])
stopifnot(names(all_dt)[1:2] == c('chrom', 'pos'))

missing_gw <- setdiff(pools_use, names(all_dt))
if (length(missing_gw) > 0) {
  stop(sprintf('pools missing from callable-all file: %s',
               paste(missing_gw, collapse = ', ')))
}

# Sum each pool column. This is the denominator that makes the reference-pool
# comparison meaningful: each pool's cluster depth is expressed relative to its
# own genome-wide callable depth, so differences in sequencing effort and in
# overall mapping rate divide out.
cat('[INFO] computing genome-wide callable depth per pool\n')
gw_depth <- sapply(pools_use, function(p) {
  sum(as.numeric(all_dt[[p]]), na.rm = TRUE)
})

# Keep the default schema exactly as it was -- step 5 reads this file, and the
# 'field' path must stay byte-identical. The extra descriptive columns are added
# only under --pool-set all, where the file is written under a different name.
gw_dt <- data.table(
  pool                   = pools_use,
  genome_callable_depth  = as.numeric(gw_depth)
)
if (pool_set == 'all') {
  gw_dt[, pool_type := ifelse(pool %in% field_pools, 'field', 'reference')]
  gw_dt[, n_callable_pos_genome := nrow(all_dt)]
  gw_dt[, mean_callable_depth := genome_callable_depth / n_callable_pos_genome]
}

cat('[INFO] genome-wide callable depth per pool:\n')
for (i in seq_len(nrow(gw_dt))) {
  cat(sprintf('         %-12s %s\n',
              gw_dt$pool[i],
              format(gw_dt$genome_callable_depth[i], big.mark = ',', scientific = FALSE)))
}


# ---- Write outputs ----------------------------------------------------------

# 'all' writes under distinct names so the GLMM inputs are never overwritten.
suffix <- if (pool_set == 'all') sprintf('_all%d', length(pools_use)) else ''
out_cds_path <- file.path(opt[['out-dir']],
                          sprintf('cds_depth_per_gene_pool%s.tsv', suffix))
out_gw_path  <- file.path(opt[['out-dir']],
                          sprintf('genome_callable_depth_per_pool%s.tsv', suffix))

fwrite(cds_agg_filt, out_cds_path, sep = '\t')
fwrite(gw_dt,        out_gw_path,  sep = '\t')

cat(sprintf('[INFO] wrote %s (%s rows)\n', out_cds_path,
            format(nrow(cds_agg_filt), big.mark = ',')))
cat(sprintf('[INFO] wrote %s (%s rows)\n', out_gw_path,
            format(nrow(gw_dt), big.mark = ',')))

if (pool_set == 'all') {
  # A preview of the quantity of interest, so a failed run is obvious from the
  # log.
  CYP6 <- c('CQUJHB009005', 'CQUJHB008986', 'CQUJHB009148',
            'CQUJHB017007', 'CQUJHB001163', 'CQUJHB007871')
  present <- intersect(CYP6, unique(cds_agg_filt$gene_id))
  cat('\n[INFO] preview: CpL1a CYP6 cluster depth relative to each pool\'s own\n')
  cat('[INFO] median gene. A ratio near 1 means no local elevation.\n')
  if (length(present) < length(CYP6)) {
    cat(sprintf('[WARN] only %d of the 6 focal paralogs survived the min-positions filter: %s\n',
                length(present), paste(present, collapse = ', ')))
  }
  if (length(present) > 0) {
    x <- copy(cds_agg_filt)
    x[, per_pos := total_depth / n_callable_pos]
    x[, ratio := per_pos / median(per_pos, na.rm = TRUE), by = pool]
    clu <- x[gene_id %in% present, .(cluster_ratio = mean(ratio)), by = pool]
    clu[, pool_type := ifelse(pool %in% field_pools, 'field', 'REFERENCE')]
    setorder(clu, pool_type, -cluster_ratio)
    for (i in seq_len(nrow(clu))) {
      cat(sprintf('         %-10s %-10s %6.3f\n',
                  clu$pool_type[i], clu$pool[i], clu$cluster_ratio[i]))
    }
    cat('[INFO] PIP/ME1 is the number to read first: it is near-conspecific with\n')
    cat('[INFO] the assembly, so a ratio near 1 says the elevation seen in the\n')
    cat('[INFO] field pools is a property of the sample and not of the reference.\n')
  }
  cat('\n[OK] step 4 complete (pool-set all). Next: copy both *_all*.tsv files to\n')
  cat('[OK] data/glmm_v4/ and run R/08_loci/09_cyp6_reference_pool_depth.R\n')
} else {
  cat('[OK] step 4 complete. Next: copy both files to data/glmm_v4/ and run\n')
  cat('[OK] R/06_depth/01_callable_depth_glmm.R\n')
}
