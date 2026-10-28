# =============================================================================
# R/07_multi_locus/14_matched_null_xtx.R
# =============================================================================
# Differentiation-matched null for the locus-level private-panel statistics.
# The width-matched null holds region size constant but not the differentiation
# on which the candidate regions were selected. Here comparison regions of each
# focal width (5,000 drawn per locus in proportion to chromosome length,
# excluding the candidate intervals and 8 Mb on either side, so that no
# comparison region is linked to a candidate) are retained when their mean
# XtXst lies within a tolerance of the focal region's, widening from +/-10% in
# 5% steps until at least 200 are retained. The matched set's chromosomes are
# reported: chr1 74-113 Mb is the only broad region of the genome as
# differentiated as the focal loci, so most matched regions come from it. This controls differentiation-driven
# ascertainment; it does NOT control ascertainment on environmental association
# itself, which only a full environment-permuted BayPass null would reach.
#
# Two statistics are tested per region, neither of which was a selection criterion:
#   (1) the correlation of Delta (panel alpha minus genome-wide panel alpha) with
#       the locus's environmental predictor
#   (2) the among-pool variance of Delta (donor-likeness)
#
# MRK in the BayPass output is the row index of allele_frequencies_v4.csv, so
# XtXst joins by row order (checked below).
#
# Manuscript: Supplementary Methods "Differentiation-matched null distribution";
#             Supplementary Results "Locus statistics against a differentiation-
#             matched null"; Results "Measuring donor ancestry and testing it
#             against a genome-wide null" and CpL1a section ("a correlation
#             exceeding that of every differentiation-matched region")
#
# Inputs:  data/site_coords_with_spatial_data_allvars.csv   (R/00_setup/01_site_predictors.R)
#          deposit: data/allele_frequencies_v4.csv          (R/01_filtering/01_allele_frequency_table.R)
#          deposit: results/baypass/v4/baypass_core_envsubset_summary_pi_xtx.out
#                                                           (cluster/04_baypass/baypass_core_envsubset.sh)
# Outputs: results/multilocus/Table_xtx_matched_null_summary.csv
#          results/multilocus/Table_xtx_matched_null_windows.csv
#
# Runtime: about 25 minutes
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")
suppressMessages(library(data.table))

AF_PATH  <- deposit_file("data/allele_frequencies_v4.csv")
XTX_PATH <- deposit_file("results/baypass/v4/baypass_core_envsubset_summary_pi_xtx.out")
need_files("data/site_coords_with_spatial_data_allvars.csv")

HI <- 0.95; LO <- 0.05
N_DRAW    <- 5000    # random windows drawn per locus before XtX matching
MIN_PANEL <- 5       # minimum panel SNPs for a window to be analysable
XTX_TOL   <- 0.10    # matching tolerance, as a fraction of the focal mean XtXst
BUFFER    <- 8e6     # comparison regions lie at least this far from every candidate
                     # interval: flank SNPs reach baseline correlation with the
                     # candidates within 8 Mb at all three loci (Figure 4)

# CpL2a is the 500 kb candidate core (13.5-14.0 Mb) in the private-panel analyses.
LOCI <- list(
  CpL1a = list(chr = CHR1, start = CPL1A_START, end = CPL1A_END,
               env = "bio05", panel = "QUI"),
  CpL2a = list(chr = CHR2, start = CPL2A_CORE_START,      end = CPL2A_CORE_END,
               env = "bio05", panel = "MOL"),
  CpL3a = list(chr = CHR3, start = CPL3A_START, end = CPL3A_END,
               env = "pis",   panel = "MOL")
)

# -----------------------------------------------------------------------------
# Load
# -----------------------------------------------------------------------------
say("reading XtX table ...")
xtx <- fread(XTX_PATH, select = c("MRK", "XtXst"), showProgress = FALSE)
say("  XtX rows: %s   MRK %s..%s", format(nrow(xtx), big.mark = ","),
    xtx$MRK[1], xtx$MRK[nrow(xtx)])

hdr   <- names(fread(AF_PATH, nrows = 0, showProgress = FALSE))
pools <- setdiff(hdr, c("chromosome", "position", "cumulative_position"))
field <- setdiff(pools, REFS)

say("reading allele frequencies ...")
af <- fread(AF_PATH, select = c("chromosome", "position", pools), showProgress = FALSE)
say("  AF rows: %s", format(nrow(af), big.mark = ","))

stopifnot(nrow(af) == nrow(xtx))
stopifnot(identical(as.numeric(xtx$MRK), as.numeric(seq_len(nrow(xtx)))))
af[, xtxst := xtx$XtXst]
rm(xtx); invisible(gc())

env <- fread("data/site_coords_with_spatial_data_allvars.csv", showProgress = FALSE)
setkey(env, site); env <- env[field]

# Panel membership depends only on the three reference frequencies, so compute
# it once genome-wide rather than per window.
q <- af[[QUI]]; p <- af[[PIP]]; m <- af[[MOL]]
ok <- is.finite(q) & is.finite(p) & is.finite(m)
PANEL <- list(
  QUI = ok & q >= HI & p <= LO & m <= LO,
  PIP = ok & p >= HI & q <= LO & m <= LO,
  MOL = ok & m >= HI & q <= LO & p <= LO
)
say("genome-wide panel sizes: QUI %s  PIP %s  MOL %s",
    format(sum(PANEL$QUI), big.mark = ","), format(sum(PANEL$PIP), big.mark = ","),
    format(sum(PANEL$MOL), big.mark = ","))

F <- as.matrix(af[, ..field])
chrom <- af$chromosome; pos <- af$position; xt <- af$xtxst
rm(af); invisible(gc())

# Per-chromosome index ranges, for fast window extraction by position
chrs  <- unique(chrom)
chrIx <- lapply(chrs, function(c) which(chrom == c))
names(chrIx) <- chrs
chrLen <- sapply(chrs, function(c) max(pos[chrIx[[c]]]))
say("chromosomes: %s", paste(sprintf("%s (%.0f Mb)", sub("ENA_OZ00", "", chrs),
    chrLen / 1e6), collapse = ", "))

# -----------------------------------------------------------------------------
# Per-window statistics
# -----------------------------------------------------------------------------
# Genome-wide background alpha per panel: panel SNPs OUTSIDE all three candidate
# intervals. The tested statistic is Delta = alpha_locus - alpha_genome, so the
# null is built on Delta too. Testing raw alpha instead inflates the null,
# because the genome-wide background is itself environmentally structured
# (alpha_genome ~ bio05 r = +0.69; ~ impervious surface r = +0.66), so EVERY window
# inherits that signal.
in_cand <- rep(FALSE, length(pos))
for (L in LOCI) in_cand <- in_cand | (chrom == L$chr & pos >= L$start & pos <= L$end)
ALPHA_GENOME <- lapply(names(PANEL), function(p) {
  rows <- which(PANEL[[p]] & !in_cand)
  if (length(rows) < MIN_PANEL) return(NULL)
  colMeans(F[rows, , drop = FALSE], na.rm = TRUE)
})
names(ALPHA_GENOME) <- names(PANEL)
say("genome background panel sizes (outside all candidate intervals): QUI %s  PIP %s  MOL %s",
    format(sum(PANEL$QUI & !in_cand), big.mark = ","),
    format(sum(PANEL$PIP & !in_cand), big.mark = ","),
    format(sum(PANEL$MOL & !in_cand), big.mark = ","))

win_stats <- function(chr, start, end, envvar, panel_name) {
  ix  <- chrIx[[chr]]
  sel <- ix[pos[ix] >= start & pos[ix] <= end]
  if (!length(sel)) return(NULL)
  pk <- sel[PANEL[[panel_name]][sel]]
  out <- list(n_snps = length(sel), n_panel = length(pk),
              mean_xtx = mean(xt[sel], na.rm = TRUE),
              r = NA_real_, r_alpha = NA_real_, slope = NA_real_, var_mol = NA_real_)
  if (length(pk) >= MIN_PANEL) {
    alpha <- colMeans(F[pk, , drop = FALSE], na.rm = TRUE)
    x     <- as.numeric(env[[envvar]])
    delta <- alpha - ALPHA_GENOME[[panel_name]]
    if (sd(delta) > 0) {
      out$r     <- cor(delta, x)          # Delta, the tested statistic
      out$slope <- coef(lm(delta ~ x))[2]
    }
    if (sd(alpha) > 0) out$r_alpha <- cor(alpha, x)   # raw alpha, for comparison
    out$var_mol <- var(delta)
  }
  out
}

is_candidate <- function(chr, start, end) {
  any(sapply(LOCI, function(L)
    L$chr == chr && start <= L$end + BUFFER && end >= L$start - BUFFER))
}

set.seed(42)   # the random comparison windows

results <- list(); summary_rows <- list()

for (nm in names(LOCI)) {
  L <- LOCI[[nm]]
  W <- L$end - L$start
  focal <- win_stats(L$chr, L$start, L$end, L$env, L$panel)
  say("")
  say("######## %s  (%.2f Mb, panel %s-private, predictor %s) ########",
      nm, W / 1e6, L$panel, L$env)
  say("  focal: %s SNPs, %d panel SNPs, mean XtXst %.2f, r = %+.3f, var(alpha) = %.3g",
      format(focal$n_snps, big.mark = ","), focal$n_panel, focal$mean_xtx,
      focal$r, focal$var_mol)

  # draw random windows of the same width, excluding the three candidate intervals
  drawn <- list(); tries <- 0
  while (length(drawn) < N_DRAW && tries < N_DRAW * 20) {
    tries <- tries + 1
    c1 <- sample(chrs, 1, prob = chrLen / sum(chrLen))
    if (chrLen[[c1]] <= W) next
    s1 <- runif(1, 1, chrLen[[c1]] - W)
    if (is_candidate(c1, s1, s1 + W)) next
    st <- win_stats(c1, s1, s1 + W, L$env, L$panel)
    if (is.null(st) || is.na(st$r)) next
    st$chr <- c1; st$start <- s1
    drawn[[length(drawn) + 1]] <- st
  }
  D <- rbindlist(lapply(drawn, as.data.table))
  say("  drew %d analysable random windows (from %d attempts)", nrow(D), tries)

  # match on differentiation
  tol <- XTX_TOL * focal$mean_xtx
  M <- D[abs(mean_xtx - focal$mean_xtx) <= tol]
  used_tol <- XTX_TOL
  while (nrow(M) < 200 && used_tol < 1) {
    used_tol <- used_tol + 0.05
    M <- D[abs(mean_xtx - focal$mean_xtx) <= used_tol * focal$mean_xtx]
  }
  say("  XtX-matched windows: %d (tolerance +/-%.0f%% of focal mean XtXst %.2f)",
      nrow(M), 100 * used_tol, focal$mean_xtx)
  say("    matched mean XtXst: median %.2f, range %.2f-%.2f",
      median(M$mean_xtx), min(M$mean_xtx), max(M$mean_xtx))
  say("    width-only (unmatched) mean XtXst: median %.2f", median(D$mean_xtx))
  comp <- M[, .N, by = .(chr = sub("_OZ.*", "", sub("ENA_OZ00431", "chr", chr)))][order(chr)]
  say("    matched regions by chromosome: %s", paste(comp$chr, comp$N, collapse = ", "))

  pct <- function(v, focal_v, side) {
    if (side == "abs") mean(abs(v) >= abs(focal_v), na.rm = TRUE)
    else if (focal_v >= 0) mean(v >= focal_v, na.rm = TRUE)
    else mean(v <= focal_v, na.rm = TRUE)
  }
  say("  --- ancestry-environment association (r) ---")
  say("    focal r = %+.3f", focal$r)
  say("    vs WIDTH-ONLY null : P(|r| >= |focal|) = %.4f   directional P = %.4f",
      pct(D$r, focal$r, "abs"), pct(D$r, focal$r, "dir"))
  say("    vs XtX-MATCHED null: P(|r| >= |focal|) = %.4f   directional P = %.4f",
      pct(M$r, focal$r, "abs"), pct(M$r, focal$r, "dir"))
  say("  --- among-pool variance in %s-likeness ---", L$panel)
  say("    focal var = %.4g ; matched median = %.4g ; ratio = %.1fx ; P = %.4f",
      focal$var_mol, median(M$var_mol, na.rm = TRUE),
      focal$var_mol / median(M$var_mol, na.rm = TRUE),
      mean(M$var_mol >= focal$var_mol, na.rm = TRUE))

  M[, locus := nm]; results[[nm]] <- M
  summary_rows[[nm]] <- data.table(
    locus = nm, width_mb = W / 1e6, panel = L$panel, env = L$env,
    focal_n_panel = focal$n_panel, focal_mean_xtx = focal$mean_xtx,
    focal_r = focal$r, focal_var = focal$var_mol,
    n_drawn = nrow(D), n_matched = nrow(M), tol_frac = used_tol,
    matched_median_xtx = median(M$mean_xtx),
    matched_chr1 = sum(M$chr == CHR1), matched_chr2 = sum(M$chr == CHR2),
    matched_chr3 = sum(M$chr == CHR3),
    widthonly_median_xtx = median(D$mean_xtx),
    p_abs_widthonly = pct(D$r, focal$r, "abs"),
    p_abs_matched   = pct(M$r, focal$r, "abs"),
    p_dir_matched   = pct(M$r, focal$r, "dir"),
    p_var_matched   = mean(M$var_mol >= focal$var_mol, na.rm = TRUE),
    var_ratio       = focal$var_mol / median(M$var_mol, na.rm = TRUE))
}

S <- rbindlist(summary_rows)
fwrite(S, "results/multilocus/Table_xtx_matched_null_summary.csv")
fwrite(rbindlist(results), "results/multilocus/Table_xtx_matched_null_windows.csv")
say("")
say("wrote results/multilocus/Table_xtx_matched_null_summary.csv")
say("wrote results/multilocus/Table_xtx_matched_null_windows.csv")
say("done.")
