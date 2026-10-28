# =============================================================================
# R/04_pop_structure/05_pool_diversity.R
# =============================================================================
# Per-pool genome-wide diversity: expected heterozygosity from pool allele
# frequencies, the mean of 2p(1 - p) over all SNPs. This is a relative index
# across pools of comparable coverage, not a coverage-corrected pi. Field-pool
# diversity is regressed on maximum summer temperature (bio05) and percent
# impervious surface; reference-pool values are printed for contrast.
#
# Manuscript: Supplementary Results "Genome-wide diversity and differentiation
#             from the reference pools"; Results (reference heterozygosities
#             PIP 0.227, MOL 0.151, QUI 0.079); He covariate of Figure S5
#
# Inputs:  data/site_coords_with_spatial_data_allvars.csv   (site predictors; no producing script in the repository)
#          data/sample_categories.csv                       (sample metadata)
#          deposit: data/allele_frequencies_v4.csv
# Outputs: results/multilocus/Table_per_pool_diversity.csv  (pool, exp_het, n_snps)
#            read by R/04_pop_structure/06_fst_vs_temperature.R,
#            R/07_multi_locus/19_background_composition.R and
#            R/07_multi_locus/22_background_ancestry_vs_environment.R
#
# Runtime: a few minutes (reads the 27-pool allele-frequency table)
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/sample_metadata.R")   # field_names

suppressMessages(library(data.table))

AF_PATH  <- deposit_file("data/allele_frequencies_v4.csv")
ENV_FILE <- "data/site_coords_with_spatial_data_allvars.csv"
need_files(ENV_FILE)

OUT <- "results/multilocus/Table_per_pool_diversity.csv"
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

env <- fread(ENV_FILE, showProgress = FALSE)

# -----------------------------------------------------------------------------
# Per-pool expected heterozygosity
# -----------------------------------------------------------------------------
hdr   <- names(fread(AF_PATH, nrows = 0, showProgress = FALSE))
pools <- setdiff(hdr, c("chromosome", "position", "cumulative_position"))
say("reading %s for genome-wide diversity ...", AF_PATH)
af <- fread(AF_PATH, select = pools, showProgress = FALSE)
say("matrix rows: %s", format(nrow(af), big.mark = ","))
M   <- as.matrix(af)
rm(af); invisible(gc())
hetv <- colMeans(2 * M * (1 - M), na.rm = TRUE)
nsnp <- colSums(!is.na(M))
rm(M); invisible(gc())
het <- data.table(pool = names(hetv), exp_het = as.numeric(hetv),
                  n_snps = as.numeric(nsnp))
fwrite(het, OUT)
say("wrote %s", OUT)

hf <- het[pool %in% field_names]
setkey(hf, pool)
hf <- hf[env[, .(site, bio05, pis)], nomatch = 0, on = c(pool = "site")]

say("")
say("=== per-pool expected heterozygosity ===")
say("field pools: %d   mean 2p(1-p) = %.4f  (range %.4f to %.4f)",
    nrow(hf), mean(hf$exp_het), min(hf$exp_het), max(hf$exp_het))
for (v in c("bio05", "pis")) {
  f <- lm(hf$exp_het ~ hf[[v]]); s <- summary(f)
  say("  het ~ %-6s slope %+.6f  adjR2 %6.3f  P = %.3g",
      v, coef(f)[2], s$adj.r.squared, coef(s)[2, 4])
}
say("")
say("=== reference pools, for contrast (colony age differs: MOL ~160 gen) ===")
print(het[pool %in% c("FL1", "ME1", "CAL1")], row.names = FALSE)

say("")
say("done.")
