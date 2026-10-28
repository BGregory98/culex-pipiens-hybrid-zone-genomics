#!/bin/bash
# =============================================================================
# cluster/08_cpl2a_pip_mol/01_cpl2a_depth_mapq.sh
# =============================================================================
# Per-base depth (samtools depth -a -q 10 -Q 20) and mapping-quality profiles
# of the three reference pools, ME1 (PIP), FL1 (QUI) and CAL1 (MOL), across
# chr2:11-16 Mb and the CpL2a interval chr2:13.5-14.0 Mb, and the runs of
# > 100 bp at zero depth in the interval, which mark mapping failures as
# distinct from conserved sequence. Ran ~15 min.
#
# Step: 11, CpL2a PIP-MOL divergence (cluster/README.md)
# Manuscript: Supplementary Methods "Indel calling and per-base depth
#             profiling at CpL2a"; Supplementary Results "PIP--MOL sequence
#             divergence at CpL2a beyond SNPs"
#
# Inputs:  samples/{ME1,FL1,CAL1}/samfile/pip_1021/<pool>.filtered.bam
#          (01_alignment/filter_index_flagstat.sh)
# Outputs (cpl2a_analysis/depth/, per pool):
#   <SAMPLE>_cpl2a_wide_depth.tsv       per-base depth, chr2:11-16 Mb
#   <SAMPLE>_cpl2a_candidate_depth.tsv  per-base depth, chr2:13.5-14.0 Mb
#   <SAMPLE>_cpl2a_mapq_dist.tsv        MAPQ distribution of reads in the interval
#   <SAMPLE>_cpl2a_flagstat.txt         flag summary for the interval
#   <SAMPLE>_cpl2a_zero_runs.tsv        zero-depth runs > 100 bp in the interval
# =============================================================================
#SBATCH --job-name=cpl2a_depth
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH -o logs/cpl2a_depth.%j.out
#SBATCH -e logs/cpl2a_depth.%j.err

set -Eeuo pipefail
trap 'echo "[ERROR] line $LINENO status $?"' ERR

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

PROJECT="${PROJECT_ROOT}"
REF="${REF_FASTA}"

# CpL2a wide region: chr2 11.0-16.0 Mb (candidate is 13.5-14.0 Mb)
# Chromosome name in BAM headers uses pipe form
CHR2_BAM="ENA|OZ004312|OZ004312.1"
REGION="${CHR2_BAM}:11000000-16000000"
REGION_NARROW="${CHR2_BAM}:13500000-14000000"

OUTDIR="${PROJECT}/cpl2a_analysis/depth"
mkdir -p "$OUTDIR" logs

source "${CONDA_SH}"
conda activate "${ENV_CPL2A}"
module load samtools 2>/dev/null || true
command -v samtools >/dev/null || { echo "[ERROR] samtools not found" >&2; exit 1; }

SAMPLES=("ME1" "FL1" "CAL1")
LABELS=("PIP" "QUI" "MOL")

for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    LABEL="${LABELS[$i]}"
    BAM="${PROJECT}/samples/${SAMPLE}/samfile/pip_1021/${SAMPLE}.filtered.bam"

    if [[ ! -f "$BAM" ]]; then
        echo "[ERROR] BAM not found: $BAM" >&2
        exit 1
    fi

    echo "=== ${SAMPLE} (${LABEL}) ==="

    # --- Per-base depth across wide region ---
    # -a includes zero-depth positions (critical for finding mapping deserts)
    # -J includes deletions in depth count
    # -q 10 minimum mapping quality (matches mpileup_v4 settings)
    # -Q 20 minimum base quality
    echo "  Extracting per-base depth (wide region)..."
    samtools depth \
        -a \
        -q 10 \
        -Q 20 \
        -r "$REGION" \
        "$BAM" \
        > "${OUTDIR}/${SAMPLE}_cpl2a_wide_depth.tsv"

    echo "  Extracting per-base depth (candidate interval)..."
    samtools depth \
        -a \
        -q 10 \
        -Q 20 \
        -r "$REGION_NARROW" \
        "$BAM" \
        > "${OUTDIR}/${SAMPLE}_cpl2a_candidate_depth.tsv"

    # --- Per-read MAPQ distribution in the candidate interval ---
    # Extract MAPQ for every read overlapping the interval
    echo "  Extracting MAPQ distribution..."
    samtools view -q 0 -F 4 "$BAM" "$REGION_NARROW" \
        | awk '{print $5}' \
        | sort -n \
        | uniq -c \
        | awk '{print $2"\t"$1}' \
        > "${OUTDIR}/${SAMPLE}_cpl2a_mapq_dist.tsv"

    # --- Flagstat for the candidate interval ---
    echo "  Running flagstat..."
    samtools view -b "$BAM" "$REGION_NARROW" \
        | samtools flagstat - \
        > "${OUTDIR}/${SAMPLE}_cpl2a_flagstat.txt"

    # --- Quick summary stats ---
    n_bases=$(wc -l < "${OUTDIR}/${SAMPLE}_cpl2a_candidate_depth.tsv")
    n_zero=$(awk '$3 == 0' "${OUTDIR}/${SAMPLE}_cpl2a_candidate_depth.tsv" | wc -l)
    mean_dp=$(awk '{s+=$3; n++} END {printf "%.1f", s/n}' "${OUTDIR}/${SAMPLE}_cpl2a_candidate_depth.tsv")
    echo "  ${LABEL}: ${n_bases} bases, ${n_zero} at zero depth, mean depth ${mean_dp}"
    echo ""
done

# --- Identify zero-depth stretches shared or private to each reference ---
echo "=== Zero-depth stretches in candidate interval ==="
for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    LABEL="${LABELS[$i]}"
    FILE="${OUTDIR}/${SAMPLE}_cpl2a_candidate_depth.tsv"

    # Find runs of consecutive zero-depth positions > 100 bp
    awk '$3 == 0 {
        if (start == 0) { start = $2; end = $2 }
        else if ($2 == end + 1) { end = $2 }
        else {
            if (end - start + 1 > 100) printf "%s\t%d\t%d\t%d\n", $1, start, end, end-start+1
            start = $2; end = $2
        }
    }
    END {
        if (start > 0 && end - start + 1 > 100) printf "%s\t%d\t%d\t%d\n", $1, start, end, end-start+1
    }' "$FILE" > "${OUTDIR}/${SAMPLE}_cpl2a_zero_runs.tsv"

    n_runs=$(wc -l < "${OUTDIR}/${SAMPLE}_cpl2a_zero_runs.tsv")
    echo "  ${LABEL} (${SAMPLE}): ${n_runs} zero-depth runs > 100 bp"
    if [[ "$n_runs" -gt 0 ]]; then
        head -10 "${OUTDIR}/${SAMPLE}_cpl2a_zero_runs.tsv" | while read -r line; do
            echo "    $line"
        done
    fi
done

echo ""
echo "[✓] Depth analysis complete. Outputs in: $OUTDIR"
