# =============================================================================
# R/helpers/ref_colors.R
# =============================================================================
#
# Reference-pool colours and display labels for plotting. CAL1 = molestus,
# FL1 = quinquefasciatus, ME1 = pipiens s.s.
#
# Source from a consumer script:
#   source('R/helpers/ref_colors.R')
#   scale_color_manual(values = ref_cols)
#   # or for specific subset orderings:
#   scale_color_manual(values = ref_cols[c('FL1', 'ME1', 'CAL1')])
# =============================================================================

ref_cols <- c(
  CAL1 = '#C27CA6',
  FL1  = '#E9C13F',
  ME1  = '#5C88E6'
)

# -----------------------------------------------------------------------------
# Display labels
# -----------------------------------------------------------------------------
# Figures label the reference pools by lineage (MOL, QUI, PIP). The labels are
# display only: the pool identifiers are data keys -- they index columns in the
# VCF, the allele-frequency matrix and every derived table.
#
# Use when rendering axis titles, facet strips, legends and annotations:
#   source('R/helpers/ref_colors.R')
#   scale_color_manual(values = ref_cols, labels = ref_labels)
#   facet_wrap(~ ref, labeller = as_labeller(ref_pair_labels))
#
# CAL1 = molestus, FL1 = quinquefasciatus, ME1 = pipiens s.s.

ref_labels <- c(
  CAL1 = 'MOL',
  FL1  = 'QUI',
  ME1  = 'PIP'
)

# Reference-pair labels for the windowed-f3 facets, in both orderings used
# across the scripts.
ref_pair_labels <- c(
  'FL1-CAL1' = 'QUI-MOL',
  'CAL1-FL1' = 'QUI-MOL',
  'ME1-CAL1' = 'PIP-MOL',
  'CAL1-ME1' = 'PIP-MOL',
  'ME1-FL1'  = 'PIP-QUI',
  'FL1-ME1'  = 'PIP-QUI'
)

# Convenience: rewrite any character vector of pool ids to display labels,
# leaving non-reference entries (field pools) untouched.
relabel_refs <- function(x) {
  x <- as.character(x)
  hit <- x %in% names(ref_labels)
  x[hit] <- unname(ref_labels[x[hit]])
  x
}
