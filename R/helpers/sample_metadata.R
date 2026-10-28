# =============================================================================
# R/helpers/sample_metadata.R
# =============================================================================
#
# Pool names. Source after config.R (working directory = repository root):
#   source("R/helpers/sample_metadata.R")
#
# Defines:
#   sample_info_all  data/sample_categories.csv (sample, sample_type, ...)
#   field_names      the 24 field pools, sorted
#   ref_names        the three reference pools: CAL1, FL1, ME1
#   all_samples      c(field_names, ref_names), the 27 pools analysed
#   QUI, PIP, MOL    reference pool for each lineage (FL1 = Cx. quinquefasciatus,
#                    ME1 = Cx. pipiens f. pipiens, CAL1 = Cx. pipiens f. molestus)
#   REFS             c(QUI, PIP, MOL)
#
# sample_categories.csv also lists pools that were sequenced but not analysed
# (FL2, FL3, ME2, ME3 and three colonies); they are not in all_samples.
# =============================================================================

sample_info_all <- read.csv("data/sample_categories.csv")

field_names <- sort(sample_info_all$sample[sample_info_all$sample_type == "field"])
ref_names   <- c("CAL1", "FL1", "ME1")
all_samples <- c(field_names, ref_names)

QUI <- "FL1"; PIP <- "ME1"; MOL <- "CAL1"; REFS <- c(QUI, PIP, MOL)
