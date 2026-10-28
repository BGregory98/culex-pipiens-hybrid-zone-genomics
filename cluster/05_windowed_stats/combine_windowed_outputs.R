# =============================================================================
# cluster/05_windowed_stats/combine_windowed_outputs.R
# =============================================================================
# Concatenates the per-task windowed f3 CSVs into one table and adds the
# one-sided lower-tail P value of each window's block-jackknife Z score, its
# Benjamini-Hochberg adjustment within each reference pairing (popB, popC;
# all 24 field populations x windows as one family) and the admixed flag
# (adjusted P < 0.05). Run once after all array tasks finish, from the
# project root: Rscript scripts/05_windowed_stats/combine_windowed_outputs.R 1000000
# The FST block combines outputs of the windowed-FST runs (archived, not used
# in the paper) and is skipped when none are present.
#
# Step: 7, windowed f3 (cluster/README.md)
# Manuscript: Methods "Windowed f3 admixture scan"; Supplementary Methods
#             "Significance scheme for the windowed f3 admixture
#             classification"
#
# Output windowed_F3_all_comparisons_v4_1M.csv is read by
# R/07_multi_locus/08_windowed_f3.R (repo data/windowed_F3_all_comparisons_v4_1M.csv).
#
# Usage:
#   Rscript combine_windowed_outputs.R [win_size_bp]
#
#   win_size_bp defaults to 1000000 (1 Mb). Combines only files matching
#   that size suffix.
#
# Inputs:
#   SequencingData/data/windowed_fst_parts_v4/windowed_FST_*_<size>.csv
#   SequencingData/data/windowed_f3_parts_v4/windowed_F3_*_<size>.csv
#
# Outputs:
#   SequencingData/data/windowed_FST_all_comparisons_v4_<size>.csv
#   SequencingData/data/windowed_F3_all_comparisons_v4_<size>.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ---- Parse arguments ------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
win_size_bp <- if (length(args) >= 1) as.integer(args[1]) else 1000000L
if (is.na(win_size_bp) || win_size_bp < 1000) {
  stop("win_size_bp must be a positive integer in bp.")
}

format_win_size <- function(n) {
  if (n %% 1e6 == 0) return(paste0(n %/% 1e6, "M"))
  if (n %% 1e3 == 0) return(paste0(n %/% 1e3, "K"))
  return(as.character(n))
}
SIZE_LABEL <- format_win_size(win_size_bp)
cat("Combining at window size:", SIZE_LABEL, "\n\n")

# ---- Combine FST outputs ---------------------------------------------------

fst_dir   <- "SequencingData/data/windowed_fst_parts_v4"
fst_glob  <- file.path(fst_dir, sprintf("windowed_FST_*_%s.csv", SIZE_LABEL))
fst_files <- Sys.glob(fst_glob)

if (length(fst_files) == 0) {
  cat("No FST files matched", fst_glob, "; skipping FST combine.\n")
} else {
  cat("FST: combining", length(fst_files), "files...\n")
  fst_all <- rbindlist(lapply(fst_files, fread), fill = TRUE)
  out_fst <- sprintf("SequencingData/data/windowed_FST_all_comparisons_v4_%s.csv",
                     SIZE_LABEL)
  fwrite(fst_all, out_fst)
  cat("  Wrote:", out_fst, "(", nrow(fst_all), "rows )\n")
  cat("  Pairs (popA, popB):\n")
  fst_all[, .N, by = .(popA, popB)] |>
    print(nrows = 5)
  cat("\n")
}

# ---- Combine F3 outputs and add BH-adjusted p-values ----------------------

f3_dir   <- "SequencingData/data/windowed_f3_parts_v4"
f3_glob  <- file.path(f3_dir, sprintf("windowed_F3_*_%s.csv", SIZE_LABEL))
f3_files <- Sys.glob(f3_glob)

if (length(f3_files) == 0) {
  cat("No F3 files matched", f3_glob, "; skipping F3 combine.\n")
} else {
  cat("F3: combining", length(f3_files), "files...\n")
  f3_all <- rbindlist(lapply(f3_files, fread), fill = TRUE)

  # Per-stratum BH adjustment: stratify by reference pair (popB, popC).
  # Within each pair, p-values across (focal, window) combinations are
  # adjusted together. This is the manuscript-relevant scope: "of all
  # (focal, window) tests for THIS admixture hypothesis, which survive?"
  f3_all[, p_one_sided := pnorm(z_score)]  # lower tail; f3* < 0 = admixture
  f3_all[, p_adjusted  := p.adjust(p_one_sided, method = "BH"),
         by = .(popB, popC)]
  f3_all[, admix_bh    := p_adjusted < 0.05]

  out_f3 <- sprintf("SequencingData/data/windowed_F3_all_comparisons_v4_%s.csv",
                    SIZE_LABEL)
  fwrite(f3_all, out_f3)
  cat("  Wrote:", out_f3, "(", nrow(f3_all), "rows )\n")
  cat("  Reference pairs (popB, popC):\n")
  f3_all[, .(n_rows = .N,
             n_admix = sum(admix_bh, na.rm = TRUE),
             pct_admix = round(100 * mean(admix_bh, na.rm = TRUE), 1)),
         by = .(popB, popC)] |>
    print()
}

cat("\n[", format(Sys.time()), "] Done.\n")
