# =============================================================================
# R/07_multi_locus/10_windowed_pool_pca.R
# =============================================================================
# Genome-wide windowed PCA of field-pool allele frequencies: the Pool-seq
# analogue of the windowed PCA with which Cebrian-Camison et al. (2026) found
# inversions 1q, 3p and 3q in individually sequenced Cx. pipiens. An inversion
# whose frequency varies among pools should appear as a contiguous run of
# windows where (a) PC1 explains an unusually large share of variance, (b) the
# pool ordering on PC1 departs from the genome-wide ordering, and (c) the
# window's PC structure is shared with its neighbours but not with the rest of
# the genome (local PCA; Li & Ralph 2019). An inversion at similar frequency in
# every pool is invisible to this approach, and karyotypes cannot be recovered.
#
# Manuscript: Figure S11 (drawn by 11_windowed_pool_pca_figure.R); Results "A
#             pre-adapted molestus inversion tracks urbanization independently
#             of temperature"; Methods "Inversion detection and breakpoint
#             confirmation"; Supplementary Methods "Windowed and local PCA
#             scan" (1,066 and 534 windows); Supplementary Results "A
#             genome-wide windowed PCA recovers the inversion without prior
#             coordinates"
#
# Inputs:  data/sample_categories.csv              (via R/helpers/sample_metadata.R)
#          data/site_coords_with_spatial_data.csv  (R/00_setup/01_site_predictors.R)
#          deposit: data/allele_counts_v4.tsv      (cluster/04_baypass/vcf_to_allele_counts.sh)
# Outputs: results/multilocus/windowed_pca/windowed_pca_1Mb.csv
#          results/multilocus/windowed_pca/windowed_pca_500kb.csv
#          results/multilocus/windowed_pca/windowed_pca_scores_1Mb.csv
#          results/multilocus/windowed_pca/windowed_pca_scores_500kb.csv
#
# Runtime: dominated by streaming the 2.0 GB allele-count table
# =============================================================================
#
# Method
# Per window (1 Mb and 500 kb, non-overlapping), SNP filters match
# compute_region_pca() in 03_haplotype_block_diagnostics.R and the null:
# depth >= 10 in all 24 field pools, field-pool mean frequency in
# [0.05, 0.95]. Each SNP is centred on its field-pool mean; the 24 x 24
# cross-product matrix of centred frequencies gives PC1 % variance identical to
# prcomp(t(centered), center = FALSE). Cross-products are additive across
# windows, so the genome-wide PCA is the eigen-decomposition of their sum. Per
# window we also accumulate the Reynolds-style per-pool FST to CAL1 (as
# compute_mean_fst_to_donor(); CAL1 depth >= 10 additionally required). The
# local-PCA distance between two windows is lostruct's pc_dist (k = 2,
# L1-normalised eigenvalues), and the windows are placed on two MDS axes with
# cmdscale().
#
# The allele-count table is headerless: chrom, pos, ref, alt, then 27 "ref:alt"
# count columns in VCF sample order (POOL_ORDER below).

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")

library(data.table)
library(tidyverse)

COUNTS_FILE <- deposit_file("data/allele_counts_v4.tsv")
ENV_FILE    <- "data/site_coords_with_spatial_data.csv"
need_files(ENV_FILE)
OUT_DIR     <- "results/multilocus/windowed_pca"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Column order of the 27 count columns (VCF sample order)
POOL_ORDER <- c("Beltsville","UFP","QP","UMBC","Hillcrest","Bowie","WCP",
                "ColPark","UM","Cst","NatHist","ElCity","ME1","FL1",
                "ColHeights","Jerald","AGCB","NorthBroad","Clarksville","SERC",
                "Cheverly","IMET","Pax","CHC","SSC","CAL1","Brighton")

MIN_DEPTH   <- 10
MAF_MIN     <- 0.05
MIN_SNPS    <- 200       # windows with fewer filtered SNPs are dropped
CHUNK_ROWS  <- 500000L
WIN_SIZES   <- c(1e6, 5e5)   # must nest: larger is a multiple of smaller
K_LOSTRUCT  <- 2L

CHR <- c(chr1 = CHR1, chr2 = CHR2, chr3 = CHR3)

field_pools <- field_names
stopifnot(length(field_pools) == 24L, all(field_pools %in% POOL_ORDER))
field_idx <- match(field_pools, POOL_ORDER)
cal_idx   <- match("CAL1", POOL_ORDER)


# =============================================================================
# Stage 1: stream the counts file, accumulate per-window summaries
# =============================================================================

# Per-window summary from filtered frequency rows of one window.
summarise_window <- function(fq, fq_cal) {
  centred <- fq - rowMeans(fq)
  fst_ok  <- !is.na(fq_cal)
  fst_sum <- vapply(seq_len(ncol(fq)), function(j) {
    p <- fq[fst_ok, j]; q <- fq_cal[fst_ok]
    pbar <- (p + q) / 2; den <- pbar * (1 - pbar)
    v <- den > 0
    sum((p[v] - q[v])^2 / den[v])
  }, numeric(1))
  fst_n <- vapply(seq_len(ncol(fq)), function(j) {
    p <- fq[fst_ok, j]; q <- fq_cal[fst_ok]
    pbar <- (p + q) / 2
    sum(pbar * (1 - pbar) > 0)
  }, numeric(1))
  list(cp = crossprod(centred), n = nrow(fq), fst_sum = fst_sum, fst_n = fst_n)
}

base_w  <- max(WIN_SIZES)
windows <- setNames(lapply(WIN_SIZES, function(w) list()), format(WIN_SIZES, scientific = FALSE))
leftover <- NULL
skip <- 0L
t0 <- Sys.time()

process_rows <- function(dt) {
  # dt: chrom, pos, freq matrix columns f1..f24, fcal
  fq     <- as.matrix(dt[, paste0("f", seq_along(field_idx)), with = FALSE])
  for (w in WIN_SIZES) {
    key <- paste(dt$chrom, floor(dt$pos / w), sep = "|")
    grp <- split(seq_len(nrow(dt)), key)
    for (g in names(grp)) {
      ii <- grp[[g]]
      if (!is.null(windows[[format(w, scientific = FALSE)]][[g]]))
        stop("window ", g, " seen twice: counts file is not position-sorted")
      windows[[format(w, scientific = FALSE)]][[g]] <<- summarise_window(fq[ii, , drop = FALSE], dt$fcal[ii])
    }
  }
}

repeat {
  chunk <- fread(COUNTS_FILE, sep = "\t", header = FALSE, skip = skip,
                 nrows = CHUNK_ROWS, colClasses = "character", showProgress = FALSE)
  n_read <- nrow(chunk)
  if (n_read == 0L) break
  skip <- skip + n_read

  chunk <- chunk[V1 %in% CHR]
  if (nrow(chunk) > 0L) {
    cnt_cols <- paste0("V", 4L + seq_along(POOL_ORDER))
    need     <- c(field_idx, cal_idx)
    alt <- dp <- matrix(NA_real_, nrow(chunk), length(need))
    for (k in seq_along(need)) {
      parts <- tstrsplit(chunk[[cnt_cols[need[k]]]], ":", fixed = TRUE)
      r <- as.numeric(parts[[1]]); a <- as.numeric(parts[[2]])
      alt[, k] <- a; dp[, k] <- r + a
    }
    fq <- alt / dp
    fq[dp < MIN_DEPTH] <- NA_real_
    nf <- length(field_idx)
    ffield <- fq[, seq_len(nf), drop = FALSE]
    keep <- rowSums(is.na(ffield)) == 0
    mf   <- rowMeans(ffield)
    keep <- keep & mf >= MAF_MIN & mf <= 1 - MAF_MIN

    dt <- data.table(chrom = chunk$V1[keep], pos = as.integer(chunk$V2[keep]))
    dt[, paste0("f", seq_len(nf)) := as.data.table(ffield[keep, , drop = FALSE])]
    dt[, fcal := fq[keep, nf + 1L]]
    dt <- rbind(leftover, dt)

    # Hold back the last (possibly incomplete) base window
    last_key <- paste(dt$chrom, floor(dt$pos / base_w))
    tail_key <- last_key[length(last_key)]
    hold     <- last_key == tail_key
    leftover <- dt[hold]
    if (any(!hold)) process_rows(dt[!hold])
  }
  cat(sprintf("[%s] %.1fM rows read, %.1f min\n", format(Sys.time(), "%H:%M:%S"),
              skip / 1e6, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (n_read < CHUNK_ROWS) break
}
if (!is.null(leftover) && nrow(leftover) > 0L) process_rows(leftover)

cache <- list(windows = windows, field_pools = field_pools,
              config = list(MIN_DEPTH = MIN_DEPTH, MAF_MIN = MAF_MIN, WIN_SIZES = WIN_SIZES))


# =============================================================================
# Stage 2: per-window statistics, genome-wide reference, local-PCA MDS
# =============================================================================

top_k <- function(cp, n, k = K_LOSTRUCT) {
  e <- eigen(cp / n, symmetric = TRUE)
  list(values = e$values[1:k], vectors = e$vectors[, 1:k, drop = FALSE],
       pct1 = 100 * e$values[1] / sum(e$values), pc1 = e$vectors[, 1])
}

# lostruct pc_dist with L1-normalised eigenvalues
pc_dist_mat <- function(eig_list) {
  m <- length(eig_list)
  lam <- lapply(eig_list, function(e) e$values / sum(abs(e$values)))
  D <- matrix(0, m, m)
  for (i in seq_len(m - 1)) for (j in (i + 1):m) {
    cross <- crossprod(eig_list[[i]]$vectors, eig_list[[j]]$vectors)^2
    d2 <- sum(lam[[i]]^2) + sum(lam[[j]]^2) - 2 * sum(outer(lam[[i]], lam[[j]]) * cross)
    D[i, j] <- D[j, i] <- sqrt(max(d2, 0))
  }
  D
}

env <- fread(ENV_FILE)[match(field_pools, site)]

analyse_size <- function(w) {
  wl  <- cache$windows[[format(w, scientific = FALSE)]]
  wl  <- wl[vapply(wl, function(x) x$n, numeric(1)) >= MIN_SNPS]
  keys <- do.call(rbind, strsplit(names(wl), "|", fixed = TRUE))
  meta <- tibble(chrom = keys[, 1], start = as.numeric(keys[, 2]) * w,
                 n_snps = vapply(wl, function(x) x$n, numeric(1))) %>%
    mutate(chr = names(CHR)[match(chrom, CHR)], mid_mb = (start + w / 2) / 1e6)
  ord  <- order(match(meta$chr, names(CHR)), meta$start)
  meta <- meta[ord, ]; wl <- wl[ord]

  gw_cp <- Reduce(`+`, lapply(wl, `[[`, "cp"))
  gw_n  <- sum(meta$n_snps)
  gw    <- top_k(gw_cp, gw_n)

  eig  <- lapply(wl, function(x) top_k(x$cp, x$n))
  fst  <- t(vapply(wl, function(x) x$fst_sum / x$fst_n, numeric(24)))

  # Orient PC1 so higher scores = more MOL-like (lower FST to CAL1)
  pc1 <- t(vapply(seq_along(eig), function(i) {
    v <- eig[[i]]$pc1
    if (cor(v, fst[i, ]) > 0) -v else v
  }, numeric(24)))

  mds <- cmdscale(pc_dist_mat(eig), k = 2)

  res <- meta %>%
    mutate(pc1_pct   = vapply(eig, `[[`, numeric(1), "pct1"),
           cor_gw    = abs(apply(pc1, 1, cor, y = gw$pc1)),
           cor_fst   = abs(vapply(seq_len(nrow(pc1)), function(i) cor(pc1[i, ], fst[i, ]), numeric(1))),
           cor_pis   = vapply(seq_len(nrow(pc1)), function(i) cor(pc1[i, ], env$pis), numeric(1)),
           mds1 = mds[, 1], mds2 = mds[, 2])

  # Orient MDS axes so CpL3a sits at positive values (sign is arbitrary)
  in3q <- res$chr == "chr3" & res$mid_mb > CPL3A_START / 1e6 & res$mid_mb < CPL3A_END / 1e6
  for (a in c("mds1", "mds2")) if (median(res[[a]][in3q]) < 0) res[[a]] <- -res[[a]]

  scores <- as_tibble(pc1, .name_repair = ~ field_pools) %>%
    bind_cols(res %>% select(chr, start, mid_mb)) %>%
    pivot_longer(all_of(field_pools), names_to = "pool", values_to = "pc1_score") %>%
    left_join(env %>% select(pool = site, pis), by = "pool")

  tag <- if (w >= 1e6) paste0(w / 1e6, "Mb") else paste0(w / 1e3, "kb")
  fwrite(res,    file.path(OUT_DIR, sprintf("windowed_pca_%s.csv", tag)))
  fwrite(scores, file.path(OUT_DIR, sprintf("windowed_pca_scores_%s.csv", tag)))
  list(res = res, scores = scores, tag = tag, gw_pct1 = gw$pct1)
}


# =============================================================================
# Stage 3: per-interval console summary
# =============================================================================

# Inversions of Cebrian-Camison et al. (2026) and this study's candidate loci.
intervals <- tribble(
  ~chr,   ~start, ~end,    ~label,       ~source,
  "chr1",  87.0,  102.0,   "1q",         "Cebrian-Camison",
  "chr3",   3.8,   18.5,   "3p",         "Cebrian-Camison",
  "chr3", 149.0,  161.0,   "3q",         "Cebrian-Camison",
  "chr1", CPL1A_START / 1e6, CPL1A_END / 1e6, "CpL1a", "This study",
  "chr2", CPL2A_START / 1e6, CPL2A_END / 1e6, "CpL2a", "This study",
  "chr3", CPL3A_START / 1e6, CPL3A_END / 1e6, "CpL3a", "This study"
)

for (w in WIN_SIZES) {
  a <- analyse_size(w)

  # Console summary per interval
  res <- a$res
  cat(sprintf("\n==== %s windows: %d windows, genome-wide PC1 = %.1f%%, PC1 95th pct = %.1f%% ====\n",
              a$tag, nrow(res), a$gw_pct1, quantile(res$pc1_pct, 0.95)))
  bg <- res
  for (i in seq_len(nrow(intervals))) {
    iv <- intervals[i, ]
    sub <- filter(res, chr == iv$chr, mid_mb > iv$start, mid_mb < iv$end)
    if (nrow(sub) == 0) next
    pct_rank <- function(x, ref) mean(ref <= median(x))
    cat(sprintf("%-6s %-16s n=%3d  PC1%%var med %.1f (genome pct %.2f)  |r gw| med %.2f  |r FST| med %.2f  r pis med %+.2f  MDS1 med %+.3f (pct %.2f)\n",
                iv$label, iv$source, nrow(sub), median(sub$pc1_pct), pct_rank(sub$pc1_pct, bg$pc1_pct),
                median(sub$cor_gw), median(sub$cor_fst), median(sub$cor_pis),
                median(sub$mds1), pct_rank(sub$mds1, bg$mds1)))
  }
}
