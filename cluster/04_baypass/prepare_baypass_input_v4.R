# =============================================================================
# cluster/04_baypass/prepare_baypass_input_v4.R
# =============================================================================
# Converts the allele-count table into BayPass inputs: interleaved REF/ALT
# counts for all 27 pools and for the 24 field pools, SNP IDs, the eight
# environmental covariates scaled and transposed (predictors as rows), and
# pool-size files (80 haploid genomes per pool).
#
# Step: 6, BayPass (cluster/README.md)
# Manuscript: Methods "Genome-environment association with BayPass"
#
# Inputs:  --counts allele_counts_v4.tsv (vcf_to_allele_counts.sh);
#          --samples data/metadata/locations_v4.txt; --categories
#          data/metadata/sample_categories.csv; --env
#          data/metadata/site_coords_with_spatial_data.csv
# Outputs (to --outdir):
#
#   - baypass_input_counts.txt          all pools, REF/ALT count pairs
#   - baypass_snp_ids.txt               SNP IDs in same order
#   - baypass_STD_counts.txt            field pools only (env-subset analysis)
#   - baypass_STD_snp_ids.txt           SNP IDs for the env subset (same as above)
#   - baypass_env_matrix.txt            scaled, transposed env covariates
#   - baypass_poolsize_allpops.txt      1 × n_pools (haploid pool sizes, all pools)
#   - baypass_poolsize_envsubset.txt    1 × n_field (haploid pool sizes, field only)
#
# Usage:
#   Rscript prepare_baypass_input_v4.R \
#     --counts     <allele_counts.tsv> \
#     --samples    <pool name list, one per line, in VCF column order> \
#     --categories <sample_categories.csv> \
#     --env        <site_coords_with_spatial_data.csv> \
#     --outdir     <directory for baypass_* output files>
#
# Published run (from the project root, R/4.5.3 module):
#   Rscript scripts/04_baypass/prepare_baypass_input_v4.R \
#     --counts     allele_counts_v4.tsv \
#     --samples    data/metadata/locations_v4.txt \
#     --categories data/metadata/sample_categories.csv \
#     --env        data/metadata/site_coords_with_spatial_data.csv \
#     --outdir     SequencingData/data/baypass_v4
# =============================================================================

suppressPackageStartupMessages({
  library(stringr)
  library(purrr)
  library(dplyr)
  library(data.table)
  library(optparse)
})

# ---- CLI args -------------------------------------------------------------

opt_list <- list(
  make_option("--counts",     type = "character", default = NULL,
              help = "Allele counts TSV from vcf_to_allele_counts.sh"),
  make_option("--samples",    type = "character", default = NULL,
              help = "Pool name list, one per line, in VCF column order"),
  make_option("--categories", type = "character", default = NULL,
              help = "sample_categories.csv (cols: sample, sample_type, sample_type_2)"),
  make_option("--env",        type = "character", default = NULL,
              help = "site_coords_with_spatial_data.csv (env covariates per field pool)"),
  make_option("--outdir",     type = "character", default = ".",
              help = "Directory for baypass_* output files [default: %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

required <- c("counts", "samples", "categories", "env")
missing  <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing) > 0) {
  stop("Missing required args: --", paste(missing, collapse = ", --"),
       "\nRun with --help for usage.")
}

ALLELE_COUNTS_PATH <- opt$counts
SAMPLES_PATH       <- opt$samples
SAMPLE_CATEGORIES  <- opt$categories
ENV_PATH           <- opt$env
OUT_DIR            <- opt$outdir
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

HAPLOID_SIZE <- 80L  # 40 individuals × 2 = 80 haploid genomes per pool

# ---- Load samples and define field subset --------------------------------

samples <- readLines(SAMPLES_PATH) |> str_trim() |> discard(~ .x == "")
cat("Samples in VCF order (n =", length(samples), "):\n  ",
    paste(samples, collapse = ", "), "\n\n")

# Canonical field-pool list. Use the same source the post-BayPass
# classifier uses to avoid drift between scripts.
sample_info <- read.csv(SAMPLE_CATEGORIES)
field_names <- sort(sample_info$sample[sample_info$sample_type == "field"])

cat("Field pools (n =", length(field_names), "):\n  ",
    paste(field_names, collapse = ", "), "\n\n")

# Sanity: every field name should be present in `samples`. Catches stale
# sample_categories.csv (e.g., still listing pools that are not in the VCF).
missing_in_vcf <- setdiff(field_names, samples)
if (length(missing_in_vcf) > 0) {
  stop("Field pools listed in ", SAMPLE_CATEGORIES, " are missing from ",
       SAMPLES_PATH, ": ",
       paste(missing_in_vcf, collapse = ", "),
       "\n  -> Update sample_categories.csv to reflect the pool design in use.")
}

# ---- Load allele counts ---------------------------------------------------

cat("Reading allele counts from:", ALLELE_COUNTS_PATH, "\n")
dat <- fread(ALLELE_COUNTS_PATH, header = FALSE, sep = "\t")

# Sanity: should be 4 metadata columns + one column per sample.
expected_ncol <- 4 + length(samples)
if (ncol(dat) != expected_ncol) {
  stop("Column-count mismatch: ", ALLELE_COUNTS_PATH, " has ", ncol(dat),
       " columns but expected ", expected_ncol,
       " (= 4 + ", length(samples), " samples).")
}

setnames(dat, 1:4, c("CHROM", "POS", "REF", "ALT"))

# Sanity: chromosome names should be in underscore form (ENA_OZ00xxxx_...).
chrom_examples <- unique(head(dat$CHROM, 200))
if (any(grepl("\\|", chrom_examples))) {
  stop("Chromosome names contain '|' — vcf_to_allele_counts.sh should have ",
       "canonicalized these to underscore form. Re-run the converter, or ",
       "check that the VCF was built against cxpip_genome.sanitized.fa.")
}
cat("Unique chromosomes seen:", paste(unique(dat$CHROM), collapse = ", "), "\n\n")

# ---- Split RD:AD into integer columns per pool ----------------------------

cat("Splitting RD:AD columns for each pool... ")
for (i in seq_along(samples)) {
  col_idx <- 4L + i
  s       <- samples[i]
  splits  <- tstrsplit(dat[[col_idx]], ":", fixed = TRUE)
  dat[, (paste0(s, "_REF")) := as.integer(splits[[1]])]
  dat[, (paste0(s, "_ALT")) := as.integer(splits[[2]])]
}

# Drop the original RD:AD string columns; keep CHROM/POS/REF/ALT and the
# 2 × n_samples integer columns.
dat[, (paste0("V", 5:(4 + length(samples)))) := NULL]
cat("done.\n\n")

# ---- Build SNP IDs --------------------------------------------------------

dat[, SNP_ID := paste(CHROM, POS, sep = "_")]

# Counts-only matrix in canonical BayPass column order: REF/ALT *interleaved*
# per pool (sample1_REF sample1_ALT sample2_REF sample2_ALT ...).
all_cols_interleaved <- unlist(lapply(samples, function(s) {
  c(paste0(s, "_REF"), paste0(s, "_ALT"))
}))
geno_counts <- dat[, ..all_cols_interleaved]

# ---- Write all-pools BayPass inputs ---------------------------------------

write.table(
  geno_counts,
  file = file.path(OUT_DIR, "baypass_input_counts.txt"),
  quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE
)
write.table(
  dat$SNP_ID,
  file = file.path(OUT_DIR, "baypass_snp_ids.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)
cat("Wrote: baypass_input_counts.txt  (",
    nrow(geno_counts), "x", ncol(geno_counts), ")\n")
cat("Wrote: baypass_snp_ids.txt        (", nrow(dat), "rows )\n\n")

# ---- Build env-subset (field pools only) ----------------------------------

subset_cols <- unlist(lapply(field_names, function(s) {
  c(paste0(s, "_REF"), paste0(s, "_ALT"))
}))

missing_cols <- setdiff(subset_cols, colnames(geno_counts))
if (length(missing_cols) > 0) {
  stop("Field-pool columns missing from geno_counts: ",
       paste(missing_cols, collapse = ", "))
}

geno_counts_env <- geno_counts[, ..subset_cols]

write.table(
  geno_counts_env,
  file = file.path(OUT_DIR, "baypass_STD_counts.txt"),
  quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE
)
write.table(
  dat$SNP_ID,
  file = file.path(OUT_DIR, "baypass_STD_snp_ids.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)
cat("Wrote: baypass_STD_counts.txt    (",
    nrow(geno_counts_env), "x", ncol(geno_counts_env), ")\n")
cat("Wrote: baypass_STD_snp_ids.txt    (same SNPs as above)\n\n")

# ---- Build environmental covariate matrix --------------------------------

env_raw <- read.csv(ENV_PATH) |>
  dplyr::select(-any_of(c("lat", "lon"))) |>
  dplyr::rename(sample = site)

# Sanity: every field pool must be in env_raw.
missing_env <- setdiff(field_names, env_raw$sample)
if (length(missing_env) > 0) {
  stop("Field pools missing from env data: ", paste(missing_env, collapse = ", "))
}

# Enforce same row order as field_names (and therefore as the count columns).
env_data <- env_raw |>
  filter(sample %in% field_names) |>
  arrange(match(sample, field_names))

env_mat <- env_data |> dplyr::select(-sample) |> as.matrix()
env_mat_scaled <- scale(env_mat)
env_mat_scaled <- t(env_mat_scaled)  # BayPass wants predictors as rows

write.table(
  env_mat_scaled,
  file = file.path(OUT_DIR, "baypass_env_matrix.txt"),
  quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE
)
cat("Wrote: baypass_env_matrix.txt    (",
    nrow(env_mat_scaled), "predictors x", ncol(env_mat_scaled), "field pools )\n")
cat("  Predictors (rows):", paste(rownames(env_mat_scaled), collapse = ", "), "\n\n")

# ---- Pool-size files ------------------------------------------------------

npop_all <- ncol(geno_counts) / 2L
write.table(
  t(rep(HAPLOID_SIZE, npop_all)),
  file = file.path(OUT_DIR, "baypass_poolsize_allpops.txt"),
  quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE
)
cat("Wrote: baypass_poolsize_allpops.txt    (",
    npop_all, "pools, all =", HAPLOID_SIZE, ")\n")

npop_env <- ncol(geno_counts_env) / 2L
write.table(
  t(rep(HAPLOID_SIZE, npop_env)),
  file = file.path(OUT_DIR, "baypass_poolsize_envsubset.txt"),
  quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE
)
cat("Wrote: baypass_poolsize_envsubset.txt  (",
    npop_env, "field pools, all =", HAPLOID_SIZE, ")\n\n")

cat("[✓] All BayPass input files written to:", OUT_DIR, "\n")
