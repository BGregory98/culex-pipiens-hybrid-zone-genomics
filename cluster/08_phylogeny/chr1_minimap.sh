#!/bin/bash
# =============================================================================
# cluster/08_phylogeny/chr1_minimap.sh
# =============================================================================
# Aligns chromosome 1 of the JHB2020 Cx. quinquefasciatus assembly
# (CM027410.1) to chromosome 1 of idCulPipi1 (OZ004311.1) with minimap2 2.30
# (-x asm10), for the CpL1a CYP6 microsynteny plot. Submitted from the
# project root.
#
# Step: 12, CYP6 phylogeny, consensus and microsynteny (cluster/README.md)
# Manuscript: Figure S9; Supplementary Methods "Coding consequences,
#             phylogenetic placement, microsynteny and consensus sequences of
#             the CpL1a CYP6 paralogs"
#
# Inputs:  reference_genomes/chr1_only/chr1_pipiens.fa, chr1_quinq.fa
#          (chromosome 1 of each assembly)
# Outputs: reference_genomes/chr1_only/chr1_quinq_vs_pip_chr1.asm10.paf (read
#          by R/08_loci/06_cyp6_microsynteny.R; deposit
#          data/reference_genome/chr1_only/chr1_quinq_vs_pip_chr1.asm10.paf)
# =============================================================================
#SBATCH --job-name=chr1_minimap
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
mkdir -p logs

MM2="software/minimap2-2.30/minimap2"

"$MM2" -x asm10 -t 4 \
  reference_genomes/chr1_only/chr1_pipiens.fa \
  reference_genomes/chr1_only/chr1_quinq.fa \
  > reference_genomes/chr1_only/chr1_quinq_vs_pip_chr1.asm10.paf
