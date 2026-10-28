# =============================================================================
# R/helpers/cluster_presets.R
# =============================================================================
#
# Configuration for the three candidate loci (cluster1/2/3): genomic interval,
# polarization references, fixation-filter threshold, and environmental
# predictor for each.
#
# Consumed by:
#   R/07_multi_locus/06_allele_frequency_clines.R
#   R/07_multi_locus/07_reference_fst_decomposition.R
#   R/07_multi_locus/08_windowed_f3.R
#
# Source from a consumer script (after R/helpers/config.R):
#   source('R/helpers/cluster_presets.R')
#   cfg <- cluster_presets[['cluster1']]   # or 'cluster2' / 'cluster3'
#
# Then pull the fields the consumer needs from cfg:
#   cluster_chr       <- cfg$chr
#   cluster_start     <- cfg$start
#   cluster_end       <- cfg$end
#   polarize_ref_A    <- cfg$polarize_A
#   polarize_ref_B    <- cfg$polarize_B
#   fix_threshold     <- cfg$fix_threshold
#   cluster_predictor <- cfg$predictor
#   cluster_pred_label <- cfg$predictor_label
#
# Notes:
#   - `fix_threshold` is the |Delta-AF(polarize_A, polarize_B)| cutoff used
#     to identify near-fixation focal SNPs. It differs across loci because
#     the AF-contrast distributions differ (chr1: enriched in fixation-class
#     SNPs, threshold can be strict at 0.95; chr2/chr3: sparser, threshold
#     relaxed to 0.70).
#   - cluster1 and cluster3 bounds are CpL1a and CpL3a from R/helpers/loci.R.
#     cluster2 is the 500 kb CpL2a candidate-SNP core (13.5-14.0 Mb), not the
#     2.5 Mb CpL2a haplotype block (12.5-15.0 Mb) defined in loci.R.
#   - Chromosome names are in underscore form (matches GFF3 + downstream
#     BayPass / AF outputs).
# =============================================================================

source("R/helpers/loci.R")
cluster_presets <- list(

  cluster1 = list(
    chr             = CHR1,
    start           = CPL1A_START,
    end             = CPL1A_END,
    polarize_A      = 'FL1',
    polarize_B      = 'CAL1',
    fix_threshold   = 0.95,
    predictor       = 'bio05',
    predictor_label = 'Maximum Summer Temperature (bio05)'
  ),

  cluster2 = list(
    chr             = CHR2,
    start           = CPL2A_CORE_START,
    end             = CPL2A_CORE_END,
    polarize_A      = 'FL1',
    polarize_B      = 'ME1',
    fix_threshold   = 0.70,
    predictor       = 'bio05',
    predictor_label = 'Maximum Summer Temperature (bio05)'
  ),

  cluster3 = list(
    chr             = CHR3,
    start           = CPL3A_START,
    end             = CPL3A_END,
    polarize_A      = 'ME1',
    polarize_B      = 'CAL1',
    fix_threshold   = 0.70,
    predictor       = 'tree_canopy',
    predictor_label = 'Percent Tree Canopy'
  )

)
