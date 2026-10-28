# =============================================================================
# R/helpers/config.R
# =============================================================================
#
# Paths for every R script. Each script sources this file first:
#
#   source("R/helpers/config.R")
#
# Two roots:
#   repository root   the working directory. Small inputs under data/ and every
#                     output under results/ and figures/ are read and written
#                     relative to it.
#   data deposit      the unpacked Dryad deposit holding the large inputs
#                     (allele counts and frequencies, VCF, BayPass output,
#                     annotation). Its location is the environment variable
#                     CULEX_DATA, set for example in a .Renviron file at the
#                     repository root:
#                       CULEX_DATA=/path/to/culex-deposit
#                     deposit_file() resolves a path inside it.
#
# A missing input stops the script with its path; nothing is skipped silently.
# =============================================================================

if (!file.exists("R/helpers/config.R")) {
  stop("Run the scripts from the repository root (the folder that holds R/ and data/), ",
       "e.g.  Rscript R/07_multi_locus/01_per_pool_ancestry.R", call. = FALSE)
}

DATA_ROOT <- Sys.getenv("CULEX_DATA", unset = "")

# Path to a file in the data deposit; stops if CULEX_DATA is unset or the file
# is not there.
deposit_file <- function(rel) {
  if (!nzchar(DATA_ROOT)) {
    stop("CULEX_DATA is not set. Point it at the unpacked data deposit, e.g. add\n",
         "  CULEX_DATA=/path/to/culex-deposit\n",
         "to a .Renviron file in the repository root and restart R.", call. = FALSE)
  }
  p <- file.path(DATA_ROOT, rel)
  if (!file.exists(p)) {
    stop(sprintf("not found in the data deposit: %s\n  (CULEX_DATA = %s)", rel, DATA_ROOT),
         call. = FALSE)
  }
  p
}

# Stops with the list of any paths that do not exist; returns them otherwise.
need_files <- function(...) {
  paths <- c(...)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("missing input(s):\n  ", paste(missing, collapse = "\n  "),
         "\n  Run the upstream script that writes it, or check the data deposit.",
         call. = FALSE)
  }
  invisible(paths)
}

# Formatted progress line: say("%d SNPs", n)
say <- function(...) cat(sprintf(...), "\n", sep = "")
