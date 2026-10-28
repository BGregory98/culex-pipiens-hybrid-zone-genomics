# =============================================================================
# R/run_all.R
# =============================================================================
# Runs the whole R stage in order, one script per R process, from the
# repository root:
#
#   Rscript R/run_all.R              # every script
#   Rscript R/run_all.R 07_multi_locus 08_loci   # only these stages
#
# Stages run in the order 00 -> 01 -> 02 -> 03 -> 04 -> 06 -> 07 -> 08 and,
# within a stage, in file-name order; each script reads only the outputs of
# scripts before it. A script's console output goes to logs/<stage>__<name>.log.
# The run stops at the first script that fails.
#
# Skipped by default (pass the stage name to include them):
#   01_filtering   rebuilds data/allele_frequencies_v4.csv from the VCF, a file
#                  the deposit already provides (about 1 h, high memory)
# Two scripts read files that are not in this study's deposit and are skipped,
# with a message, when those files are absent; their outputs are in the
# repository, so later scripts still run:
#   R/00_setup/01_site_predictors.R           the CHELSA climatologies
#                                             (docs/DATA_DEPOSIT.md)
#   R/08_loci/11_chemosensory_gene_reference.R  the curated chemosensory
#                                             annotation of Menna et al. (2026)
#
# Runtime: about 3 h in total on a 16 GB laptop; the depth GLMM
# (06_depth/01) takes about 1 h, pairwise FST (04/01) and the two matched
# nulls (07/14, 07/15) 15-30 min each.
# =============================================================================

source("R/helpers/config.R")

STAGES  <- c("00_setup", "01_filtering", "02_qc", "03_baypass", "04_pop_structure",
             "06_depth", "07_multi_locus", "08_loci")
DEFAULT_SKIP <- "01_filtering"

args   <- commandArgs(trailingOnly = TRUE)
stages <- if (length(args)) intersect(STAGES, args) else setdiff(STAGES, DEFAULT_SKIP)
if (!length(stages)) stop("no known stage among: ", paste(args, collapse = ", "), call. = FALSE)

in_deposit <- function(rel) nzchar(DATA_ROOT) && file.exists(file.path(DATA_ROOT, rel))
OPTIONAL <- list(
  "R/00_setup/01_site_predictors.R" =
    c("data/chelsa_layers/CHELSA_bio05_1981-2010_V.2.1.tif", "CHELSA rasters"),
  "R/08_loci/11_chemosensory_gene_reference.R" =
    c("data/menna_on_pip.gff3_polished", "Menna et al. (2026) annotation")
)

dir.create("logs", showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")

for (st in stages) {
  for (f in sort(list.files(file.path("R", st), pattern = "\\.R$", full.names = TRUE))) {
    if (f %in% names(OPTIONAL) && !in_deposit(OPTIONAL[[f]][1])) {
      say("skip  %s (%s not in the deposit; outputs are in the repository)", f, OPTIONAL[[f]][2])
      next
    }
    log <- file.path("logs", sprintf("%s__%s.log", st, sub("\\.R$", "", basename(f))))
    t0  <- Sys.time()
    rc  <- system2(rscript, f, stdout = log, stderr = log)
    say("%-4s  %-55s %6.0f s", if (rc == 0) "ok" else "FAIL", f,
        as.numeric(difftime(Sys.time(), t0, units = "secs")))
    if (rc != 0) stop("failed: ", f, " (see ", log, ")", call. = FALSE)
  }
}
say("done")
