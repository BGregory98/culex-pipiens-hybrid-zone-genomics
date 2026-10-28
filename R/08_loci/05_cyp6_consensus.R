# =============================================================================
# R/08_loci/05_cyp6_consensus.R
# =============================================================================
# Tests whether the consensus sequences of the six CpL1a CYP6 paralogs shift
# toward a reference lineage along the temperature gradient. Pairwise percent
# identity between each field pool's consensus and each reference pool's
# consensus is modeled with a linear mixed model of identity on centered bio05,
# reference and their interaction, with random intercepts for gene and pool;
# the interaction is tested by likelihood-ratio comparison with the additive
# model. Run at the nucleotide (whole-gene) and amino-acid (translated CDS)
# level, with per-gene, per-reference OLS slopes for description.
#
# Manuscript: Supplementary Results "The CpL1a CYP6 cluster: phylogenetic
#             placement, microsynteny and sequence similarity" (nucleotide and
#             amino-acid models); Results "Cx. quinquefasciatus ancestry and
#             sequencing read depth at a CYP6 cluster rise with summer
#             temperature" (consensus sentence); Supplementary Methods "Coding
#             consequences, phylogenetic placement, microsynteny and consensus
#             sequences of the CpL1a CYP6 paralogs"
#
# Inputs:  data/consensus/p450s/pairwise_gene_nt_identity.tsv   (cluster/08_phylogeny/p450_consensus/)
#          data/consensus/p450s/pairwise_aa_identity.tsv        (as above)
#          data/site_coords_with_spatial_data.csv               (R/00_setup/01_site_predictors.R)
# Outputs: none (model summaries printed to the console)
#
# Runtime: < 1 min
# =============================================================================

source("R/helpers/config.R")
source("R/helpers/loci.R")
source("R/helpers/sample_metadata.R")

library(data.table)  # fread()
library(dplyr)
library(lme4)        # lmer()
library(lmerTest)    # p-values / ANOVA for lmer() models

SITE_FILE   <- "data/site_coords_with_spatial_data.csv"
IDENTITY_FILES <- c(
  nucleotide = "data/consensus/p450s/pairwise_gene_nt_identity.tsv",
  amino_acid = "data/consensus/p450s/pairwise_aa_identity.tsv"
)
need_files(SITE_FILE, IDENTITY_FILES)

# -----------------------------
# bio05 at each field site
# -----------------------------

# bio05 as the mean within a 500 m buffer around each site (Table S1).
env_data <- read.csv(SITE_FILE)
env_data <- env_data[env_data$site %in% field_names, ]
env_data <- env_data %>% dplyr::rename(sample = site)

# Field pools ordered by bio05 (maximum temperature of the warmest month); used
# as the factor levels for sample.
bio05_order <- env_data %>% arrange(bio05) %>% pull(sample)

# -----------------------------
# Identity models, one per sequence level
# -----------------------------

fit_identity <- function(path, level) {
  cat(sprintf("\n==================== %s identity ====================\n", level))

  data <- fread(path) |> as.data.frame()
  data$gene = as.factor(data$gene)
  data$seq1 = as.factor(data$seq1)
  data$seq2 = as.factor(data$seq2)

  # Keep only reference-vs-field comparisons (drop field-vs-field and
  # ref-vs-ref pairs from the all-vs-all table).
  data_ref_field <- data %>%
    filter(
      (seq1 %in% ref_names & seq2 %in% field_names) |
      (seq2 %in% ref_names & seq1 %in% field_names)
    ) %>%
    mutate(
      ref_sample = ifelse(seq1 %in% ref_names, as.character(seq1), as.character(seq2)),
      sample = ifelse(seq1 %in% field_names, as.character(seq1), as.character(seq2))
    )

  data_ref_field = data_ref_field %>%
    dplyr::select(gene, ref_sample, sample, compared_sites, pct_identity)

  # Genes in chromosome order; field pools ordered along the bio05 gradient.
  data_ref_field$gene <- factor(data_ref_field$gene, levels = CYP6_GENES)
  data_ref_field$sample = factor(data_ref_field$sample, levels = bio05_order)

  df_seq <- data_ref_field %>%
    left_join(
      env_data %>% dplyr::select(sample, bio05),
      by = "sample"
    ) %>%
    mutate(
      gene = factor(gene, levels = CYP6_GENES),
      # CAL1 (MOL) is the baseline reference level.
      ref_sample = factor(ref_sample, levels = c("CAL1", "FL1", "ME1")),
      sample = factor(sample),
      # Mean-centered bio05, in degrees C: the intercept is the expected
      # identity at average bio05.
      bio05_c = as.numeric(scale(bio05, center = TRUE, scale = FALSE))
    )

  str(df_seq)

  # Percent identity on centered bio05, reference and their interaction, with
  # random intercepts for gene and sample. The interaction tests whether the
  # temperature slope differs by reference.
  model_seq <- lmer(
    pct_identity ~ bio05_c * ref_sample + (1 | gene) + (1 | sample),
    data = df_seq
  )

  # Reduced model: one temperature slope for all references.
  model_seq_no_interaction <- lmer(
    pct_identity ~ bio05_c + ref_sample + (1 | gene) + (1 | sample),
    data = df_seq
  )

  # Likelihood-ratio test of the interaction.
  print(anova(model_seq_no_interaction, model_seq))
  print(summary(model_seq))
  print(anova(model_seq))

  # Per-gene, per-reference OLS slopes of identity on bio05 (descriptive).
  gene_ref_slopes <- df_seq %>%
    group_by(gene, ref_sample) %>%
    group_modify(~{
      fit <- lm(pct_identity ~ bio05_c, data = .x)
      tibble(
        slope = coef(fit)[["bio05_c"]],
        p_val = summary(fit)$coefficients["bio05_c", "Pr(>|t|)"],
        r_squared = summary(fit)$adj.r.squared
      )
    }) %>%
    ungroup()

  print(gene_ref_slopes, n = Inf)
  invisible(model_seq)
}

for (level in names(IDENTITY_FILES)) fit_identity(IDENTITY_FILES[[level]], level)
