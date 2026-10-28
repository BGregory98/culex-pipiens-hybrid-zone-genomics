#!/bin/bash
# =============================================================================
# cluster/08_phylogeny/run_p450_phylo.sh
# =============================================================================
# Maximum-likelihood phylogeny of the CpL1a CYP6 proteins. Each protein set is
# aligned two ways with MAFFT 7.526, L-INS-i (--localpair) and G-INS-i
# (--globalpair), each also trimmed with trimAl -automated1, and IQ-TREE 3.1.1
# (-st AA -m MFP -alrt 1000 -B 1000) is run on all four alignments.
# Protein sets:
#   focal12v2  the 6 paralogs from JHB2020 and their 6 lifted copies
#   all23      focal12 + 11 reference P450s
#   all26      focal12 + 14 reference P450s: the tree in Figure S8 (L-INS-i;
#              the other three all26 variants give the same topology; model
#              chosen by BIC: LG+I+G4 for all four)
#
# Step: 12, CYP6 phylogeny, consensus and microsynteny (cluster/README.md)
# Manuscript: Figure S8; Supplementary Methods "Coding consequences,
#             phylogenetic placement, microsynteny and consensus sequences of
#             the CpL1a CYP6 paralogs"
#
# Inputs (under p450_cluster6/):
#   03_protein/focal12.prot.fa        translated CDS of the 6 paralogs, JHB2020 + lifted
#   07_refs/reduced_clean_refs.fa     curated UniProt reference P450s (11 when all23
#                                     was built, 14 when all26 was built)
#   03_protein/all26.shortnames.prot.fa  all26.prot.fa with short sequence names
#                                     (e.g. JHB_007871, B0WFB8_Cxq_CYP71B38); the
#                                     renaming command was not saved
# Outputs (under p450_cluster6/):
#   04_align/<set>.{linsi,ginsi}[.automated1].fa
#   05_trees/<set>.{linsi,ginsi}[.automated1].{treefile,contree,iqtree,log,...}
#   (repo copies: data/consensus/p450s/trees/*.treefile)
#
# The IQ-TREE runs set no -seed; the seeds IQ-TREE drew are in each .log
# (all26: linsi 545798, linsi.automated1 336885, ginsi 68010,
# ginsi.automated1 420177). Add -seed <value> to repeat a tree exactly.
# =============================================================================
#SBATCH --job-name=p450_phylo
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=16:00:00
#SBATCH --output=06_logs/p450_phylo.%j.out
#SBATCH --error=06_logs/p450_phylo.%j.err

set -euo pipefail

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

cd "${PROJECT_ROOT}/p450_cluster6"

# Load software; adjust module names to match your cluster
module purge
module load mafft
module load trimal
module load iqtree
module load seqkit

# The published alignments and trees used 4 threads (MAFFT "4 thread(s)",
# IQ-TREE "-nt 4" in their logs).
THREADS=4

mkdir -p 04_align 05_trees 06_logs

# Build combined dataset: 12 focal + 11 refs = 23 proteins
# (reduced_clean_refs.fa held 11 proteins when all23 was built; it was then
# extended to the 14 of all26, so this line now gives 26.)
cat 03_protein/focal12.prot.fa 07_refs/reduced_clean_refs.fa > 03_protein/all23.prot.fa

echo "== Sequence counts =="
echo "focal12:"
grep -c '^>' 03_protein/focal12.prot.fa
echo "all23:"
grep -c '^>' 03_protein/all23.prot.fa

echo "== Duplicate headers in all23 =="
seqkit seq -n 03_protein/all23.prot.fa | sort | uniq -d || true

run_phylo () {
    STEM="$1"
    INFA="$2"

    echo "=== Working on $STEM ==="

    mafft --thread "$THREADS" --localpair --maxiterate 1000 "$INFA" \
        > "04_align/${STEM}.linsi.fa" \
        2> "06_logs/${STEM}.linsi.mafft.log"

    mafft --thread "$THREADS" --globalpair --maxiterate 1000 "$INFA" \
        > "04_align/${STEM}.ginsi.fa" \
        2> "06_logs/${STEM}.ginsi.mafft.log"

    trimal -in "04_align/${STEM}.linsi.fa" \
           -out "04_align/${STEM}.linsi.automated1.fa" \
           -automated1 \
           > "06_logs/${STEM}.linsi.automated1.trimal.log" 2>&1

    trimal -in "04_align/${STEM}.ginsi.fa" \
           -out "04_align/${STEM}.ginsi.automated1.fa" \
           -automated1 \
           > "06_logs/${STEM}.ginsi.automated1.trimal.log" 2>&1

    for ALN in \
        "04_align/${STEM}.linsi.fa" \
        "04_align/${STEM}.linsi.automated1.fa" \
        "04_align/${STEM}.ginsi.fa" \
        "04_align/${STEM}.ginsi.automated1.fa"
    do
        PREFIX="05_trees/$(basename "${ALN%.fa}")"

        iqtree \
            -s "$ALN" \
            -st AA \
            -m MFP \
            -alrt 1000 \
            -B 1000 \
            -nt "$THREADS" \
            -pre "$PREFIX" \
            -redo
    done
}

# focal12v2 is a rerun of an earlier focal12 set under a new name, so the
# earlier focal12 files stay untouched
run_phylo focal12v2 03_protein/focal12.prot.fa
run_phylo all23     03_protein/all23.prot.fa

# all26 (Figure S8), run after reduced_clean_refs.fa was extended to 14
# reference proteins. The same concatenation as all23; the alignments and
# trees were made from the copy with shortened names.
cat 03_protein/focal12.prot.fa 07_refs/reduced_clean_refs.fa > 03_protein/all26.prot.fa
echo "all26:"
grep -c '^>' 03_protein/all26.prot.fa
run_phylo all26     03_protein/all26.shortnames.prot.fa

echo "Done."
