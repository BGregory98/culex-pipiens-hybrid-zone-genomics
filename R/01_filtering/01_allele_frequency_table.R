# =============================================================================
# R/01_filtering/01_allele_frequency_table.R
# =============================================================================
# Builds the per-pool allele-frequency table that most downstream scripts read.
# Reads the biallelic, no-missing PoolSNP VCF and takes each pool's ALT allele
# frequency from FORMAT/FREQ (PoolSNP's AD / DP, two decimals). Adds a
# cumulative genome coordinate: position plus the summed lengths of the
# preceding chromosomes. One row per SNP in VCF order, pool columns in VCF
# sample order.
#
# Manuscript: Methods "Read processing, alignment, and variant calling" (the
#             11,629,890 biallelic SNPs with complete allele-frequency data)
#
# Inputs:  deposit: data/vcfs/v4/poolsnp.biallelic.nomiss.vcf.gz
#                   (cluster/_toplevel/PoolSNP_v4.sh, then
#                    bcftools view -m2 -M2 -v snps)
# Outputs: data/allele_frequencies_v4.csv   (the deposit holds the same file
#          as data/allele_frequencies_v4.csv; downstream scripts read that copy)
#
# The companion count table, deposit: data/allele_counts_v4.tsv (no header;
# chromosome, position, REF, ALT, then one RD:AD column per pool in the same
# pool order), is written from the same VCF on the cluster by
# cluster/04_baypass/vcf_to_allele_counts.sh.
#
# Runtime: tens of minutes; holds the whole VCF FREQ matrix in memory (several GB).
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")

suppressPackageStartupMessages({
  library(data.table)
  library(VariantAnnotation)
})

VCF_PATH <- deposit_file("data/vcfs/v4/poolsnp.biallelic.nomiss.vcf.gz")
OUT_CSV  <- "data/allele_frequencies_v4.csv"

# Per-chromosome offsets for cumulative_position: the lengths of chr1 and of
# chr1 + chr2 in idCulPipi1.
OFFSETS <- setNames(c(0, 128522754, 128522754 + 213114244), c(CHR1, CHR2, CHR3))

# ---- Read VCF and extract FREQ ----------------------------------------------

say("[%s] Reading VCF: %s", format(Sys.time()), VCF_PATH)
param <- ScanVcfParam(geno = "FREQ", fixed = NA, info = NA)
vcf   <- readVcf(VCF_PATH, "csq", param = param)

# FREQ matrix: rows = SNPs, columns = pools.
freq_mat <- geno(vcf)$FREQ
storage.mode(freq_mat) <- "numeric"

say("  SNPs:  %d", nrow(freq_mat))
say("  Pools: %d", ncol(freq_mat))
say("  Pool names: %s", paste(colnames(freq_mat), collapse = ", "))

# ---- Build output table -----------------------------------------------------

gr    <- rowRanges(vcf)
chrom <- as.character(seqnames(gr))
pos   <- start(gr)

# Chromosome names in underscore form (a no-op for this VCF, whose names are
# already ENA_OZ00xxxx_OZ00xxxx.1).
chrom <- gsub("\\|", "_", chrom)

missing_chr <- setdiff(unique(chrom), names(OFFSETS))
if (length(missing_chr) > 0) {
  stop("Chromosomes in VCF missing from OFFSETS: ",
       paste(missing_chr, collapse = ", "))
}

cumulative_position <- pos + OFFSETS[chrom]

dat <- data.table(
  chromosome          = chrom,
  position            = pos,
  cumulative_position = cumulative_position
)
# Per-pool FREQ columns, in VCF pool order.
for (p in colnames(freq_mat)) {
  dat[, (p) := freq_mat[, p]]
}

# ---- Checks -----------------------------------------------------------------

n_total   <- nrow(dat)
n_with_na <- sum(rowSums(is.na(freq_mat)) > 0)
say("  Total SNPs in VCF:     %d", n_total)
say("  SNPs with any NA FREQ: %d (0 expected: the VCF has no missing pools)", n_with_na)

# ---- Write ------------------------------------------------------------------

# CRLF line ends and R's default number formatting (scipen 0, so e.g. a
# cumulative position of 1,000,000 is written 1e+06) reproduce the deposited
# table byte for byte.
fwrite(dat, OUT_CSV, eol = "\r\n", scipen = 0L)
say("Wrote %s (%d SNPs x %d pools)", OUT_CSV, n_total, ncol(freq_mat))
