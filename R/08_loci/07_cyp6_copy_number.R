# =============================================================================
# R/08_loci/07_cyp6_copy_number.R
# =============================================================================
# Implied copy number at the CpL1a CYP6 cluster and the per-paralog
# depth-temperature slopes. Depth per callable position at each gene is divided
# by the median of that quantity across all genes in the same pool (1.0 = a
# typical gene in that pool); the cluster ratio pools the six paralogs. The
# ratio is regressed on bio05 across the 24 field pools, for the cluster and
# for each paralog, and converted to copy number on a diploid baseline.
#
# Manuscript: Supplementary Results "Implied copy number and per-paralog depth
#             slopes"
#
# Inputs:  data/glmm_v4/cds_depth_per_gene_pool.tsv            (cluster/06_depth/glmm/04_aggregate_cds_depth.R)
#          data/site_coords_with_spatial_data_allvars.csv      (tracked input; no producing script in R/)
#          results/chr1/Table_chr1_cds_snp_attribution.csv     (02_candidate_snp_attribution.R)
# Outputs: results/chr1/Table_cyp6_per_paralog_depth.csv
#          results/chr1/Table_cyp6_cluster_depth_ratio.csv
#
# Runtime: < 1 min
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")

suppressMessages(library(data.table))

DEPTH_FILE <- "data/glmm_v4/cds_depth_per_gene_pool.tsv"
ENV_FILE   <- "data/site_coords_with_spatial_data_allvars.csv"
ATT_FILE   <- "results/chr1/Table_chr1_cds_snp_attribution.csv"
need_files(DEPTH_FILE, ENV_FILE, ATT_FILE)

d   <- fread(DEPTH_FILE, showProgress = FALSE)
env <- fread(ENV_FILE, showProgress = FALSE)
att <- fread(ATT_FILE, showProgress = FALSE)

# The six paralogs, in the order of the attribution table (gene-body density).
CYP6 <- att[gene_id %in% CYP6_GENES, gene_id]
stopifnot(setequal(CYP6, CYP6_GENES))
say("CYP6 paralogs: %d  (%s)", length(CYP6), paste(CYP6, collapse = ", "))
say("depth table: %s rows, %d genes, %d pools",
    format(nrow(d), big.mark = ","), uniqueN(d$gene_id), uniqueN(d$pool))

# depth per callable position, per gene per pool
d <- d[n_callable_pos > 0]
d[, dpp := total_depth / n_callable_pos]

# within-pool background: median depth-per-position across all genes
bg <- d[, .(bg_dpp = median(dpp, na.rm = TRUE)), by = pool]
d  <- merge(d, bg, by = "pool")
d[, ratio := dpp / bg_dpp]      # 1.0 = a typical gene in that pool

pools <- intersect(unique(d$pool), env$site)
setkey(env, site); e <- env[pools]
x <- as.numeric(e$bio05)
say("pools matched to bio05: %d", length(pools))

# --- cluster-level ratio per pool (pooled across the six paralogs) ---
cl <- d[gene_id %in% CYP6 & pool %in% pools,
        .(dpp = sum(total_depth) / sum(n_callable_pos)), by = pool]
cl <- merge(cl, bg, by = "pool")
cl[, ratio := dpp / bg_dpp]
setkey(cl, pool); cl <- cl[pools]

cool <- order(x)[1:6]; warm <- order(x, decreasing = TRUE)[1:6]
say("")
say("=== cluster depth ratio (1.0 = typical gene in the same pool) ===")
say("  6 coolest pools : mean %.3f   implied copy number %.2f", mean(cl$ratio[cool]), 2*mean(cl$ratio[cool]))
say("  6 warmest pools : mean %.3f   implied copy number %.2f", mean(cl$ratio[warm]), 2*mean(cl$ratio[warm]))
say("  all 24 pools    : range %.3f to %.3f", min(cl$ratio), max(cl$ratio))
f <- lm(cl$ratio ~ x); s <- summary(f)
say("  ratio ~ bio05   : slope %+.4f per C, adjR2 %.3f, P = %.3g",
    coef(f)[2], s$adj.r.squared, s$coefficients[2, 4])
say("  fitted ratio    : %.3f at coolest site -> %.3f at warmest",
    predict(f, data.frame(x = min(x))), predict(f, data.frame(x = max(x))))
say("  implied copy no.: %.2f -> %.2f  (diploid baseline of 2 assumed)",
    2*predict(f, data.frame(x = min(x))), 2*predict(f, data.frame(x = max(x))))
say("  relative rise across the gradient: %+.1f%%",
    100*(predict(f, data.frame(x = max(x)))/predict(f, data.frame(x = min(x))) - 1))

# --- per-paralog slopes ---
say("")
say("=== per-paralog depth ratio against bio05 ===")
say("%-14s %8s %8s %10s %9s %11s", "gene", "cool", "warm", "slope/C", "adjR2", "P")
res <- list()
for (g in CYP6) {
  dd <- d[gene_id == g & pool %in% pools]
  setkey(dd, pool); dd <- dd[pools]
  if (any(!is.finite(dd$ratio))) next
  f2 <- lm(dd$ratio ~ x); s2 <- summary(f2)
  say("%-14s %8.3f %8.3f %+10.4f %9.3f %11.3g", g,
      mean(dd$ratio[cool]), mean(dd$ratio[warm]),
      coef(f2)[2], s2$adj.r.squared, s2$coefficients[2, 4])
  res[[g]] <- data.table(gene_id = g, cool = mean(dd$ratio[cool]),
                         warm = mean(dd$ratio[warm]), slope = coef(f2)[2],
                         adj_r2 = s2$adj.r.squared, p = s2$coefficients[2, 4])
}
R <- rbindlist(res)
say("")
say("  paralogs with a positive slope        : %d of %d", sum(R$slope > 0), nrow(R))
say("  positive and significant at P < 0.05  : %d", sum(R$slope > 0 & R$p < 0.05))
fwrite(R, "results/chr1/Table_cyp6_per_paralog_depth.csv")
fwrite(cl[, .(pool, bio05 = x, ratio, implied_cn = 2*ratio)],
       "results/chr1/Table_cyp6_cluster_depth_ratio.csv")
say("")
say("wrote results/chr1/Table_cyp6_per_paralog_depth.csv")
say("wrote results/chr1/Table_cyp6_cluster_depth_ratio.csv")
