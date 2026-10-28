# =============================================================================
# R/08_loci/14_kdr_frequency.R
# =============================================================================
# Per-pool frequency of kdr L1014F (Vgsc, chromosome 2 position 191,086,281)
# from the raw allele counts, joined to site coordinates, environmental
# predictors and Figure 1 site numbers. Tests the frequency against bio05 and
# percent impervious surface (OLS), between metropolitan areas (Wilcoxon), and
# against genome-wide donor ancestry and geography (Spearman, with
# leave-one-out ranges for the ancestry correlations).
#
# Manuscript: Figure S7 (data, drawn by 17_kdr_figure.R); Results "A genome
#             scan identifies three candidate regions with distinct
#             architectures" (mean 0.34, range 0.08-0.82; P = 0.32 and 0.16);
#             Supplementary Results "A characterized resistance allele outside
#             the candidate set" (rho = -0.52 with PIP-private ancestry)
#
# Inputs:  results/resistance/vgsc_coding_variant_counts.tsv       (13_resistance_target_sites.R)
#          data/site_coords_with_spatial_data.csv                  (R/00_setup/01_site_predictors.R)
#          results/multilocus/per_pool_ancestry_per_pool.csv       (R/07_multi_locus/01_per_pool_ancestry.R)
# Outputs: results/resistance/Table_kdr_L1014F_by_site.csv
#
# Runtime: < 1 min
# =============================================================================
#
# The genome scan found no association of L1014F with either environmental
# axis. That is a null against the chosen predictors, not evidence of no
# spatial structure, so the frequency is also tested against metropolitan area,
# donor ancestry and geography (longitude, elevation).
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")
source("R/helpers/geo.R")

suppressPackageStartupMessages({
  library(data.table)
})

COUNTS_TSV <- 'results/resistance/vgsc_coding_variant_counts.tsv'
ENV_FILE   <- 'data/site_coords_with_spatial_data.csv'
ANC_FILE   <- 'results/multilocus/per_pool_ancestry_per_pool.csv'
need_files(COUNTS_TSV, ENV_FILE, ANC_FILE)

KDR_POS    <- 191086281L                    # L1014F, canonical CQUJHB005371.R8259
MIN_DEPTH  <- 10L                           # same threshold as 13_resistance_target_sites.R

# ---- (1) Per-pool frequency at L1014F ---------------------------------------

cnt   <- fread(COUNTS_TSV)
pools <- setdiff(names(cnt), c('chr', 'pos', 'ref', 'alt'))
row   <- cnt[chr == CHR2 & pos == KDR_POS]
stopifnot(nrow(row) == 1L)

parse_ad <- function(x) {
  s  <- tstrsplit(x, ':', fixed = TRUE)
  rd <- as.numeric(s[[1]]); ad <- as.numeric(s[[2]])
  dp <- rd + ad
  list(freq = fifelse(dp >= MIN_DEPTH, ad / dp, NA_real_), dp = dp, alt = ad)
}

per_pool <- rbindlist(lapply(pools, function(p) {
  v <- parse_ad(row[[p]])
  data.table(pool = p, alt_reads = v$alt, depth = v$dp, kdr = v$freq)
}))

refs <- per_pool[pool %in% ref_names]
cat('\n=== kdr L1014F, reference pools ===\n')
print(refs[, .(pool, depth, kdr = round(kdr, 3))])

# ---- (2) Join field pools to coordinates ------------------------------------

env <- fread(ENV_FILE)
dat <- merge(per_pool[pool %in% field_names], env,
             by.x = 'pool', by.y = 'site', all.x = TRUE)
setorder(dat, pool)
stopifnot(nrow(dat) == length(field_names), !anyNA(dat$lat))

# Site numbers as printed in manuscript Figure 1: the vector in
# R/00_setup/04_figure1_site_map.R is aligned to site_coords sorted by name and
# filtered to the field pools, which is exactly the ordering of `dat` here.
dat[, site_num := c(24,  9,  8, 16,  2,  5, 17,  1, 10,  4, 18, 19,
                    21, 15,  3, 23, 11, 13,  7, 12, 22,  6, 20, 14)]

fwrite(dat[, .(site = pool, site_num, lat, lon, alt_reads, depth,
               kdr_L1014F = round(kdr, 4), pis = round(pis, 2),
               bio05 = round(bio05, 3))],
       'results/resistance/Table_kdr_L1014F_by_site.csv')

cat(sprintf('\nField pools: n = %d, mean %.3f, range %.3f-%.3f, median depth %.0f\n',
            nrow(dat), mean(dat$kdr), min(dat$kdr), max(dat$kdr), median(dat$depth)))

for (v in c('bio05', 'pis')) {
  m <- summary(lm(dat$kdr ~ dat[[v]]))
  cat(sprintf('  kdr ~ %-6s slope %+.4f  P = %.3f  adj R2 = %.3f\n',
              v, coef(m)[2, 1], coef(m)[2, 4], m$adj.r.squared))
}
# Metro area: nearer of the two city centres, the rule used for the
# locus-level metropolitan-area terms (R/07_multi_locus/18_within_metro_models.R).
dat[, metro := fifelse(hav(lat, lon, BAL[1], BAL[2]) < hav(lat, lon, DC[1], DC[2]),
                       'Baltimore', 'Washington')]
cat(sprintf('  Baltimore  n = %2d  mean %.3f\n  Washington n = %2d  mean %.3f\n  Wilcoxon P = %.3f\n',
            sum(dat$metro == 'Baltimore'),  mean(dat[metro == 'Baltimore']$kdr),
            sum(dat$metro == 'Washington'), mean(dat[metro == 'Washington']$kdr),
            wilcox.test(kdr ~ metro, data = dat)$p.value))

# ---- (3) Ancestry and geography ----------------------------------------------
#
# L1014F is near-fixed in the molestus reference (0.91) and effectively absent
# from the other two, so the obvious candidate is not an environmental variable
# at all but each pool's molestus ancestry. alpha_genome in the per-pool ancestry
# table is that pool's mean private-panel frequency outside the focal locus, so
# it is the genome-wide background for that donor. Longitude and elevation are
# also tested, because the map reads as a southeast-high / northwest-low axis.

anc <- fread(ANC_FILE)
bg  <- dcast(anc[panel %in% c('CAL1_private', 'FL1_private', 'ME1_private')],
             pool ~ panel, value.var = 'alpha_genome')
setnames(bg, c('CAL1_private', 'FL1_private', 'ME1_private'),
             c('alpha_MOL', 'alpha_QUI', 'alpha_PIP'), skip_absent = TRUE)
dat <- merge(dat, bg, by = 'pool', all.x = TRUE)

cat('\n=== kdr against ancestry and geography ===\n')
for (v in c('alpha_MOL', 'alpha_QUI', 'alpha_PIP', 'lon', 'elevation')) {
  ok <- is.finite(dat[[v]])
  ct <- cor.test(dat$kdr[ok], dat[[v]][ok], method = 'spearman', exact = FALSE)
  m  <- summary(lm(dat$kdr[ok] ~ dat[[v]][ok]))
  cat(sprintf('  kdr ~ %-10s  Spearman rho %+.3f  P = %.4f   adj R2 = %.3f\n',
              v, ct$estimate, ct$p.value, m$adj.r.squared))
}

# Leave-one-out, because n = 24 and these are post hoc: the PIP relationship
# holds throughout (P 0.004-0.028), the MOL one does not (P 0.012-0.126).
for (v in c('alpha_PIP', 'alpha_MOL')) {
  ps <- vapply(seq_len(nrow(dat)), function(i)
          cor.test(dat$kdr[-i], dat[[v]][-i], method = 'spearman',
                   exact = FALSE)$p.value, numeric(1))
  cat(sprintf('  leave-one-out %-10s P range %.4f-%.4f\n', v, min(ps), max(ps)))
}

cat('\nWrote results/resistance/Table_kdr_L1014F_by_site.csv\n')
