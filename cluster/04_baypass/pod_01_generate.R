#!/usr/bin/env Rscript
# =============================================================================
# cluster/04_baypass/pod_01_generate.R
# =============================================================================
# Simulates a pseudo-observed dataset (POD) of 100,000 neutral SNPs for
# calibrating XtXst, following Gautier (2015) and the BayPass manual
# ("Calibration of the XtX statistic"), with simulate.baypass() from
# BayPass's baypass_utils.R. Seeded (--seed, default 42).
#
# Step: 6, BayPass (cluster/README.md)
# Manuscript: Methods "Genome-environment association with BayPass"
#
# Inputs:  --omega, --beta_params (baypass_core_envsubset.sh); --counts,
#          --poolsize (prepare_baypass_input_v4.R); --baypass_utils
#          (BAYPASS_UTILS in cluster/config.sh)
# Outputs: <outdir>/G.btapods, poolsize.btapods, pod_generation.log
#
# Strategy:
#   - Use the Ω matrix and (a_π, b_π) parameters estimated from the core-model
#     run as the neutral model.
#   - Simulate NSNP_POD pool-allele counts under that model, using the
#     observed read-coverage distribution as the binomial sample size.
#   - Run BayPass core on the POD (separate SLURM script) to get a properly
#     calibrated XtX null distribution.
#
# Usage:
#   Rscript pod_01_generate.R \
#     --omega         <core_prefix>_mat_omega.out \
#     --beta_params   <core_prefix>_summary_beta_params.out \
#     --counts        baypass_STD_counts.txt \
#     --poolsize      baypass_poolsize_envsubset.txt \
#     --baypass_utils <path>/baypass_public-master/utils/baypass_utils.R \
#     --outdir        <directory for POD files> \
#     [--nsnp 100000] [--seed 42] [--pi-maf 0.05]
#
# Published run (from the project root, R/4.5.3 module):
#   Rscript scripts/04_baypass/pod_01_generate.R \
#     --omega         SequencingData/data/baypass_v4/baypass_core_envsubset_mat_omega.out \
#     --beta_params   SequencingData/data/baypass_v4/baypass_core_envsubset_summary_beta_params.out \
#     --counts        SequencingData/data/baypass_v4/baypass_STD_counts.txt \
#     --poolsize      SequencingData/data/baypass_v4/baypass_poolsize_envsubset.txt \
#     --baypass_utils software/baypass_public-master/utils/baypass_utils.R \
#     --outdir        SequencingData/data/baypass_v4/pod
#
# Outputs (to --outdir):
#   G.btapods                simulated count file (BayPass STD format)
#   poolsize.btapods         pool size file (copy of input)
#   pod_generation.log       record of choices and dimensions
# =============================================================================

suppressPackageStartupMessages({
  library(stats)
  library(optparse)
})

# ---- CLI args -------------------------------------------------------------

opt_list <- list(
  make_option("--omega",         type = "character", default = NULL,
              help = "Ω matrix output from core run (..._mat_omega.out)"),
  make_option("--beta_params",   type = "character", default = NULL,
              help = "β_π params output from core run (..._summary_beta_params.out)"),
  make_option("--counts",        type = "character", default = NULL,
              help = "BayPass-formatted counts file (e.g., baypass_STD_counts.txt)"),
  make_option("--poolsize",      type = "character", default = NULL,
              help = "Pool-size file (e.g., baypass_poolsize_envsubset.txt)"),
  make_option("--baypass_utils", type = "character", default = NULL,
              help = "Path to baypass_utils.R (from BayPass distribution)"),
  make_option("--outdir",        type = "character", default = NULL,
              help = "Directory for POD output files"),
  make_option("--nsnp",          type = "integer",   default = 100000L,
              help = "Number of SNPs to simulate [default: %default]"),
  make_option("--seed",          type = "integer",   default = 42L,
              help = "RNG seed [default: %default]"),
  make_option("--pi_maf",        type = "double",    default = 0.05,
              help = "MAF threshold for simulate.baypass [default: %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# optparse converts --beta_params → opt$beta_params, etc.
required   <- c("omega", "beta_params", "counts", "poolsize", "baypass_utils", "outdir")
required_r <- gsub("-", "_", required)
missing    <- required[vapply(required_r, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing) > 0) {
  stop("Missing required args: --", paste(missing, collapse = ", --"),
       "\nRun with --help for usage.")
}

OMEGA_PATH    <- opt$omega
BETA_PATH     <- opt$beta_params
COUNTS_PATH   <- opt$counts
POOLSIZE_PATH <- opt$poolsize
BAYPASS_UTILS <- opt$baypass_utils
OUTDIR        <- opt$outdir
NSNP_POD      <- opt$nsnp
SEED          <- opt$seed
PI_MAF        <- opt$pi_maf

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

logfile <- file(file.path(OUTDIR, "pod_generation.log"), open = "wt")
sink(logfile, type = "output", split = TRUE)

cat("=== POD generation ===\n")
cat("Started:", format(Sys.time()), "\n")
cat("NSNP_POD =", NSNP_POD, "  SEED =", SEED, "  PI_MAF =", PI_MAF, "\n\n")

# ---- Load BayPass utilities (provides geno2YN and simulate.baypass) ----

if (!file.exists(BAYPASS_UTILS)) {
  stop("Cannot find BayPass utilities at: ", BAYPASS_UTILS)
}
source(BAYPASS_UTILS)
cat("Sourced BayPass utilities from", BAYPASS_UTILS, "\n")

# ---- Read Ω matrix ----

omega <- as.matrix(read.table(OMEGA_PATH))
npop  <- nrow(omega)
cat("Loaded Ω from", OMEGA_PATH, "  dim:", paste(dim(omega), collapse = " x "), "\n")
cat("Inferred number of pools:", npop, "\n")

# ---- Read β_π parameters ----

beta.tbl     <- read.table(BETA_PATH, header = TRUE)
pi.beta.coef <- beta.tbl$Mean[1:2]
cat("Loaded β_π parameters: a_pi =", round(pi.beta.coef[1], 4),
    "  b_pi =", round(pi.beta.coef[2], 4), "\n")

# ---- Load original counts to get NN (read-coverage matrix) ----

cat("\nLoading original count file from", COUNTS_PATH, "...\n")
orig.data <- geno2YN(COUNTS_PATH)
NN <- orig.data$NN
cat("Original NN dim:", paste(dim(NN), collapse = " x "), "\n")

# simulate.baypass() expects sample.size as nsnp x npop.
if (ncol(NN) == npop) {
  cat("NN orientation: rows = SNPs, cols = pops (correct for simulate.baypass).\n")
} else if (nrow(NN) == npop) {
  cat("NN orientation: rows = pops, cols = SNPs — transposing.\n")
  NN <- t(NN)
} else {
  stop("NN dimensions (", paste(dim(NN), collapse = " x "),
       ") do not match Ω pool count (", npop, ").")
}
total.snps <- nrow(NN)
cat("Total SNPs in original data:", total.snps, "\n")

# ---- Subsample NN rows to NSNP_POD ----

if (total.snps < NSNP_POD) {
  stop("--nsnp (", NSNP_POD, ") exceeds total SNP count (", total.snps, ").")
}
set.seed(SEED)
sample.idx <- sample.int(total.snps, NSNP_POD)
NN.sub <- NN[sample.idx, , drop = FALSE]
cat("Subsampled NN to", nrow(NN.sub), "SNPs.\n")
cat("Coverage summary across subsampled SNPs (median per pool):\n")
print(round(apply(NN.sub, 2, median), 1))

# ---- Simulate POD ----
#
# simulate.baypass() writes G.<suffix> to the *current working directory*,
# so chdir to OUTDIR for the call and restore afterward.

cat("\nCalling simulate.baypass()...\n")
old.wd <- getwd()
setwd(OUTDIR)
sim <- simulate.baypass(
  omega.mat         = omega,
  nsnp              = NSNP_POD,
  sample.size       = NN.sub,
  beta.pi           = pi.beta.coef,
  pi.maf            = PI_MAF,
  suffix            = "btapods",
  remove.fixed.loci = TRUE
)
setwd(old.wd)
cat("Simulation complete.\n")
cat("Output count file:", file.path(OUTDIR, "G.btapods"), "\n")

# ---- Copy pool size file (POD uses the same pool sizes as the real data) ----

poolsize.dst <- file.path(OUTDIR, "poolsize.btapods")
if (file.copy(POOLSIZE_PATH, poolsize.dst, overwrite = TRUE)) {
  cat("Copied pool size file to", poolsize.dst, "\n")
} else {
  warning("Failed to copy pool size file from ", POOLSIZE_PATH)
}

cat("\nFinished:", format(Sys.time()), "\n")
sink(type = "output")
close(logfile)
