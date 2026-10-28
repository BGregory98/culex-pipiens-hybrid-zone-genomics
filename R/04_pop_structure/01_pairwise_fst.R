# =============================================================================
# R/04_pop_structure/01_pairwise_fst.R
# =============================================================================
# Genome-wide pairwise FST among all 27 pools (24 field + 3 reference) from the
# PoolSNP VCF, with poolfstat's analysis-of-variance estimator and each pool
# treated as 80 chromosomes. Writes the full pairwise table and its 24 field x
# 3 reference slice; downstream scripts read these rather than re-running
# vcf2pooldata.
#
# Manuscript: Table S2; Results "mean pairwise FST 0.004, range 0.0003-0.0110"
#             (Figure S3 legend); input to Figures S3, S5 and S6
#
# Inputs:  data/sample_categories.csv                   (sample metadata)
#          deposit: data/vcfs/v4/poolsnp.biallelic.nomiss.vcf.gz
# Outputs: data/pairwise_fst_all.csv                     (pop1, pop2, fst; 351 pairs)
#          data/genome_wide_FST_for_field_vs_references.csv (sample, ref, gw_fst; 72 rows)
#
# Runtime: tens of minutes (vcf2pooldata over 11.6 M SNPs)
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/sample_metadata.R")   # field_names, ref_names

suppressPackageStartupMessages({
  library(poolfstat)
  library(VariantAnnotation)  # scanVcfHeader
  library(dplyr)
  library(tidyr)
  library(tibble)
})

VCF_FILE          <- deposit_file("data/vcfs/v4/poolsnp.biallelic.nomiss.vcf.gz")
OUT_PAIRWISE_ALL  <- "data/pairwise_fst_all.csv"
OUT_GW_REFS       <- "data/genome_wide_FST_for_field_vs_references.csv"

POOL_SIZE_HAPLOID <- 80L         # 40 diploid individuals per pool

# ---- 1. Pool order --------------------------------------------------------
# Pool names in VCF header order: vcf2pooldata assigns them to the count
# columns by position, so a sorted list would mislabel the pools.
vcf_pool_names <- samples(scanVcfHeader(VCF_FILE))

# Every reference and field pool must be in the VCF.
missing_refs   <- setdiff(ref_names,   vcf_pool_names)
missing_fields <- setdiff(field_names, vcf_pool_names)
if (length(missing_refs) > 0) {
  stop("Reference pools missing from VCF: ", paste(missing_refs, collapse = ", "))
}
if (length(missing_fields) > 0) {
  stop("Field pools missing from VCF: ", paste(missing_fields, collapse = ", "))
}

# ---- 2. Build pooldata -----------------------------------------------------
cat("[", format(Sys.time()), "] Reading VCF + building pooldata...\n")
pooldata <- vcf2pooldata(
  vcf.file  = VCF_FILE,
  poolsizes = rep(POOL_SIZE_HAPLOID, length(vcf_pool_names)),
  poolnames = vcf_pool_names
)

# ---- 3. Pairwise FST across all 27 pools ----------------------------------
cat("[", format(Sys.time()), "] Computing pairwise FST across all pools...\n")
pairwise_fst        <- compute.pairwiseFST(pooldata)
pairwise_fst_matrix <- pairwise_fst@PairwiseFSTmatrix

# ---- 4. Full pairwise table (long format) ---------------------------------
# Use the lower triangle so each unordered pair appears once.
idx <- which(lower.tri(pairwise_fst_matrix), arr.ind = TRUE)
pairwise_long <- tibble(
  pop1 = rownames(pairwise_fst_matrix)[idx[, 1]],
  pop2 = colnames(pairwise_fst_matrix)[idx[, 2]],
  fst  = pairwise_fst_matrix[idx]
)
stopifnot(nrow(pairwise_long) == choose(nrow(pairwise_fst_matrix), 2))

cat("\n[ sanity: full pairwise FST ]\n")
cat("  Pools:                ", nrow(pairwise_fst_matrix), "\n")
cat("  Unique pairs (n=C(n,2)):", nrow(pairwise_long), "\n")
cat("  FST range:            ", sprintf("%.4f - %.4f\n",
                                        min(pairwise_long$fst),
                                        max(pairwise_long$fst)))

# ---- 5. Slice the 24 field x 3 reference block ----------------------------
gw_block <- pairwise_fst_matrix[field_names, ref_names]

gw_long <- as.data.frame(gw_block) %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(cols = -sample, names_to = "ref", values_to = "gw_fst") %>%
  arrange(ref, sample)

cat("\n[ sanity: field x reference slice ]\n")
cat("  Field pools (n =", length(field_names), "):", paste(field_names, collapse = ", "), "\n\n")
cat("  Reference pools (n =", length(ref_names), "):", paste(ref_names, collapse = ", "), "\n\n")
cat("  Output rows (should be 24 x 3 = 72):", nrow(gw_long), "\n")
cat("  FST range per reference:\n")
gw_long %>%
  group_by(ref) %>%
  summarise(min = min(gw_fst), median = median(gw_fst), max = max(gw_fst)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ---- 6. Write -------------------------------------------------------------
cat("\n[", format(Sys.time()), "] Writing:", OUT_PAIRWISE_ALL, "\n")
write.csv(pairwise_long, OUT_PAIRWISE_ALL, row.names = FALSE)

cat("[", format(Sys.time()), "] Writing:", OUT_GW_REFS, "\n")
write.csv(gw_long, OUT_GW_REFS, row.names = FALSE)

# ---- 7. Field-vs-field summary statistics ---------------------------------
field_pairs <- pairwise_long %>%
  dplyr::filter(pop1 %in% field_names & pop2 %in% field_names)

cat("\n[ field-vs-field summary ]\n")
cat(sprintf('  n_pairs = %d, mean = %.4f, range = %.4f - %.4f\n',
            nrow(field_pairs), mean(field_pairs$fst),
            min(field_pairs$fst), max(field_pairs$fst)))
cat('  Min pair:  '); print(field_pairs %>% arrange(fst) %>% head(1),         row.names = FALSE)
cat('  Max pair:  '); print(field_pairs %>% arrange(desc(fst)) %>% head(1),   row.names = FALSE)

cat("[", format(Sys.time()), "] Done.\n")
