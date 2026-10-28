# =============================================================================
# R/08_loci/15_kdr_cpl1a_matched_null.R
# =============================================================================
# Tests the kdr-CpL1a association against empirical nulls matched so that the
# comparison regions are as extreme as CpL1a on whatever CpL1a was ascertained
# for. A bare Spearman correlation cannot be reported for CpL1a, because CpL1a
# was chosen for extreme differentiation and extreme environmental association,
# and both inflate its correlation with anything that tracks population
# structure -- which kdr does (14_kdr_frequency.R).
#
# Two arms, because the two halves of the CpL1a signal need different nulls.
#
#   ARM A -- ancestry. The differentiation-matched window set of
#   R/07_multi_locus/14_matched_null_xtx.R: random 2.5 Mb windows matched to
#   CpL1a on mean XtXst and panel size. For each, QUI-private alpha per pool is
#   recomputed and correlated with kdr.
#
#   ARM B -- copy number. The depth analysis uses the genome-wide depth
#   relationship as its own null, so the analogue here is: correlate every
#   gene's per-pool depth ratio with kdr and ask where the CYP6 cluster falls,
#   against contiguous six-gene runs genome-wide, and against the most variable
#   1% of runs, because the CYP6 cluster is unusually variable and a more
#   variable run has more room to correlate with anything.
#
# Both arms report the raw correlation and the partial that controls bio05 and
# genome-wide pipiens background, computed identically for focal and null, so the
# comparison is like for like.
#
# Manuscript: Supplementary Results "A characterized resistance allele outside
#             the candidate set" (13,553 runs; rho = +0.47, P = 0.022; partial
#             r = +0.56, P = 0.006; most-variable 1%, P = 0.007; 2.5 Mb window
#             P = 0.21); Supplementary Methods "Targeted analysis of
#             insecticide-resistance target sites"
#
# Inputs:  results/multilocus/Table_xtx_matched_null_windows.csv   (R/07_multi_locus/14_matched_null_xtx.R)
#          data/glmm_v4/cds_depth_per_gene_pool_all27.tsv           (cluster/06_depth/glmm/04b_aggregate_cds_depth_all27.slurm)
#          results/resistance/Table_kdr_L1014F_by_site.csv          (14_kdr_frequency.R)
#          results/multilocus/per_pool_ancestry_per_pool.csv        (R/07_multi_locus/01_per_pool_ancestry.R)
#          results/chr1/Table_chr1_gene_intersection.csv            (01_cpl1a_gene_density.R)
#          deposit: data/allele_frequencies_v4.csv                  (QUI-private panel, genome-wide)
#          deposit: data/quinq_on_pip.gff3_polished                 (gene coordinates for the contiguous-run null)
# Outputs: results/chr1/Table_kdr_cpl1a_matched_null.csv
#
# Runtime: ~5 min (frequency table)
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")

suppressPackageStartupMessages({
  library(data.table)
})

WIN_FILE   <- 'results/multilocus/Table_xtx_matched_null_windows.csv'
DEPTH_FILE <- 'data/glmm_v4/cds_depth_per_gene_pool_all27.tsv'
KDR_FILE   <- 'results/resistance/Table_kdr_L1014F_by_site.csv'
ANC_FILE   <- 'results/multilocus/per_pool_ancestry_per_pool.csv'
GENES_FILE <- 'results/chr1/Table_chr1_gene_intersection.csv'
need_files(WIN_FILE, DEPTH_FILE, KDR_FILE, ANC_FILE, GENES_FILE)
AF_PATH  <- deposit_file('data/allele_frequencies_v4.csv')
GFF_PATH <- deposit_file('data/quinq_on_pip.gff3_polished')

HI  <- 0.95; LO <- 0.05
CPL1A <- LOCI$CpL1a
set.seed(20260731)

# ---- covariates, held fixed across focal and null ---------------------------

kdrT <- fread(KDR_FILE)
anc  <- fread(ANC_FILE)
gwP  <- anc[panel == 'ME1_private', .(pool, gwPIP = alpha_genome)]
cov  <- merge(kdrT[, .(pool = site, kdr = kdr_L1014F, bio05)], gwP, by = 'pool')
setorder(cov, pool)
cov  <- cov[pool %in% field_names]
POOLS <- cov$pool
stopifnot(length(POOLS) == 24L)
Z <- as.matrix(cov[, .(bio05, gwPIP)])

# r and partial-r of an arbitrary per-pool vector against kdr
score <- function(x) {
  if (!all(is.finite(x)) || sd(x) == 0) return(c(r = NA_real_, pr = NA_real_))
  c(r  = unname(cor(cov$kdr, x, method = 'spearman')),
    pr = unname(cor(resid(lm(cov$kdr ~ Z)), resid(lm(x ~ Z)))))
}
emp_p <- function(focal, null) {              # two-sided empirical p, +1 correction
  null <- null[is.finite(null)]
  (sum(abs(null) >= abs(focal)) + 1) / (length(null) + 1)
}

# =============================================================================
# ARM A -- ancestry, against the differentiation-matched window set
# =============================================================================

# The matched windows lie on all three chromosomes, so the QUI-private panel is
# read genome-wide and each window is looked up on its own chromosome.
armA <- NULL
{
  af <- fread(AF_PATH, select = c('chromosome', 'position', QUI, PIP, MOL, POOLS),
              showProgress = FALSE)

  # QUI-private panel: fixed in the quinquefasciatus reference, absent from both others
  panel <- af[get(QUI) >= HI & get(PIP) <= LO & get(MOL) <= LO]
  cat(sprintf('ARM A: %d SNPs genome-wide, %d QUI-private panel SNPs\n', nrow(af), nrow(panel)))
  rm(af); invisible(gc())

  alpha_in <- function(chr, lo, hi) {
    p <- panel[chromosome == chr & position >= lo & position < hi]
    if (nrow(p) < 5L) return(NULL)
    vapply(POOLS, function(q) mean(p[[q]], na.rm = TRUE), numeric(1))
  }

  focalA <- score(alpha_in(CPL1A$chr, CPL1A$start, CPL1A$end))
  win <- fread(WIN_FILE)[locus == 'CpL1a']
  cat(sprintf('ARM A: %d differentiation-matched null windows\n', nrow(win)))
  nullA <- rbindlist(lapply(seq_len(nrow(win)), function(i) {
    a <- alpha_in(win$chr[i], win$start[i], win$start[i] + 2.5e6)
    if (is.null(a)) return(NULL)
    as.data.table(as.list(score(a)))
  }))
  armA <- data.table(
    arm = 'ancestry (alpha_QUI)', null_set = 'XtX + panel-size matched, 2.5 Mb',
    n_null = nrow(nullA),
    focal_r = round(focalA['r'], 3),  p_r  = round(emp_p(focalA['r'],  nullA$r), 4),
    focal_pr = round(focalA['pr'], 3), p_pr = round(emp_p(focalA['pr'], nullA$pr), 4))
}

# =============================================================================
# ARM B -- CYP6 copy number, against a variance-matched gene-level null
# =============================================================================

d <- fread(DEPTH_FILE)
d <- d[pool %in% POOLS & n_callable_pos > 0]
# Depth per callable position over the MEDIAN gene of the same pool. The median
# denominator is deliberate: a depth-weighted mean gives materially different
# slopes and would contradict the published per-paralog values
# (see 08_cyp6_depth_checks.R).
d[, per_pos := total_depth / n_callable_pos]
d[, ratio   := per_pos / median(per_pos, na.rm = TRUE), by = pool]

w <- dcast(d, gene_id ~ pool, value.var = 'ratio')
w <- w[complete.cases(w)]
setcolorder(w, c('gene_id', POOLS))
M <- as.matrix(w[, ..POOLS])
rownames(M) <- w$gene_id
cat(sprintf('ARM B: %d genes x %d pools\n', nrow(M), ncol(M)))

# exclude the whole CpL1a interval from the null, not just the six paralogs
gi <- fread(GENES_FILE)
interval_genes <- unique(gi[[grep('gene', names(gi), value = TRUE)[1]]])
cat(sprintf('ARM B: excluding %d genes in the CpL1a interval from the null\n',
            sum(rownames(M) %in% interval_genes)))

sc <- t(apply(M, 1, score))
is_focal <- rownames(M) %in% CYP6_GENES
in_locus <- rownames(M) %in% interval_genes

cluster_vec <- colMeans(M[is_focal, , drop = FALSE])
focalB <- score(cluster_vec)
focal_var <- var(cluster_vec)

# The null must be CONTIGUOUS six-gene runs, not random scattered sixes. The CYP6
# paralogs are a ~28 kb tandem array, so their depths covary and the cluster mean
# keeps that variance; averaging six genes from across the genome cancels it, and
# variance-matched scattered sixes essentially do not exist (none in 20,000
# draws). Contiguous runs reproduce the tandem-array structure by construction.
suppressPackageStartupMessages(library(rtracklayer))
g   <- import.gff3(GFF_PATH)
g   <- g[g$type %in% c('gene', 'protein_coding_gene')]
gc_ <- unique(data.table(gene_id = sub('^gene[:-]', '', as.character(g$ID)),
                         chr = as.character(seqnames(g)), start = start(g)), by = 'gene_id')
gc_ <- gc_[gene_id %in% rownames(M)]
setorder(gc_, chr, start)
ord   <- gc_$gene_id
Mo    <- M[ord, , drop = FALSE]
chrv  <- gc_$chr
runs <- lapply(seq_len(nrow(Mo) - 5L), function(i) {
  idx <- i:(i + 5L)
  if (length(unique(chrv[idx])) > 1L) return(NULL)      # no runs across chromosomes
  if (any(ord[idx] %in% interval_genes)) return(NULL)   # exclude the whole CpL1a interval
  idx
})
runs <- Filter(Negate(is.null), runs)
nullB <- rbindlist(lapply(runs, function(idx) {
  v <- colMeans(Mo[idx, , drop = FALSE])
  as.data.table(as.list(c(score(v), vr = var(v))))
}))
nullB <- nullB[is.finite(r)]
cat(sprintf('ARM B: %d contiguous 6-gene runs outside CpL1a\n', nrow(nullB)))
cat(sprintf('ARM B: CYP6 cluster among-pool variance %.4f, null median %.4f (percentile %.1f)\n',
            focal_var, median(nullB$vr), 100 * mean(nullB$vr < focal_var)))

armB <- data.table(
  arm = 'copy number (CYP6 cluster)', null_set = 'contiguous 6-gene runs, genome-wide',
  n_null = nrow(nullB),
  focal_r = round(focalB['r'], 3),  p_r  = round(emp_p(focalB['r'],  nullB$r), 4),
  focal_pr = round(focalB['pr'], 3), p_pr = round(emp_p(focalB['pr'], nullB$pr), 4))

# Variance is not itself a confound for a correlation -- correlation is
# scale-free -- but a run with little among-pool signal has its correlation
# attenuated by measurement noise, so the CYP6 cluster's exceptional variance
# does work in its favour. NO run in the genome is as variable as it is (see the
# percentile printed above), so an exact variance match does not exist. The
# closest achievable control is the most variable 1% of runs.
cut1 <- quantile(nullB$vr, 0.99, na.rm = TRUE)
nullBv <- nullB[vr >= cut1]
armBv <- data.table(
  arm = '  same, vs most-variable 1% of runs',
  null_set = sprintf('runs with variance >= %.4f (no run reaches %.4f)', cut1, focal_var),
  n_null = nrow(nullBv),
  focal_r = round(focalB['r'], 3),  p_r  = round(emp_p(focalB['r'],  nullBv$r), 4),
  focal_pr = round(focalB['pr'], 3), p_pr = round(emp_p(focalB['pr'], nullBv$pr), 4))

# Per-paralog, against all single genes outside the interval
single_null <- sc[!in_locus, ]
armB_single <- rbindlist(lapply(CYP6_GENES, function(g) {
  if (!g %in% rownames(sc)) return(NULL)
  data.table(arm = paste0('  ', g), null_set = 'all genes outside CpL1a',
             n_null = sum(is.finite(single_null[, 'r'])),
             focal_r = round(sc[g, 'r'], 3),
             p_r  = round(emp_p(sc[g, 'r'],  single_null[, 'r']), 4),
             focal_pr = round(sc[g, 'pr'], 3),
             p_pr = round(emp_p(sc[g, 'pr'], single_null[, 'pr']), 4))
}))

res <- rbindlist(list(armA, armB, armBv, armB_single), fill = TRUE)
cat('\n=== kdr-CpL1a against matched nulls ===\n')
print(res, row.names = FALSE)
fwrite(res, 'results/chr1/Table_kdr_cpl1a_matched_null.csv')

cat('\nWrote results/chr1/Table_kdr_cpl1a_matched_null.csv\n')
