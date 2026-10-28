#!/bin/bash
# =============================================================================
# cluster/08_cpl2a_pip_mol/02_cpl2a_indels_sv.sh
# =============================================================================
# Indel calling and structural-variant extraction for ME1 (PIP), CAL1 (MOL),
# and FL1 (QUI) at the CpL2a region on chromosome 2. Ran ~30 min.
#
# Step: 11, CpL2a PIP-MOL divergence (cluster/README.md)
# Manuscript: Supplementary Methods "Indel calling and per-base depth
#             profiling at CpL2a"; Supplementary Results "PIP--MOL sequence
#             divergence at CpL2a beyond SNPs"
#
# Inputs:  samples/{ME1,FL1,CAL1}/samfile/pip_1021/<pool>.filtered.bam
#          (01_alignment/filter_index_flagstat.sh); reference: REF_FASTA;
#          discordant_pair_analysis/delly/<pool>_inversions.bcf
#          (07_inversion/run_delly_all3.sh)
#
# Part A: bcftools indel calling
#   Runs mpileup + call on the three reference BAMs at chr2:11-16 Mb,
#   retaining both SNPs and indels. Filters to indels only, then extracts
#   the CpL2a candidate interval (13.5-14.0 Mb).
#
# Part B: DELLY chr2 extraction
#   DELLY was already run genome-wide on all three references (inversions).
#   This script extracts chr2 calls from the existing BCF output without
#   re-running DELLY.
#
# Part C: Discordant read pairs at CpL2a
#   Extracts the paired primary reads in the candidate interval from each BAM
#   and counts those with the mate on another chromosome or an insert
#   > 10 kb, signatures of structural variants.
#
# Outputs (cpl2a_analysis/):
#   indels/cpl2a_indels_all.vcf.gz          all variants (SNPs and indels), chr2:11-16 Mb
#   indels/cpl2a_indels_candidate.tsv       candidate-interval indels, tabular
#   indels/cpl2a_indels_pip_mol_fixed.tsv   indels fixed-different PIP vs MOL
#   delly/<SAMPLE>_inversions_chr2_all.tsv, <SAMPLE>_inversions_cpl2a.tsv
#                                           DELLY calls on chr2 / at 11-16 Mb
#   discordant/<SAMPLE>_cpl2a_discordant.bam  paired primary reads in the interval
# =============================================================================
#SBATCH --job-name=cpl2a_indel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -o logs/cpl2a_indel.%j.out
#SBATCH -e logs/cpl2a_indel.%j.err

set -Eeuo pipefail
trap 'echo "[ERROR] line $LINENO status $?"' ERR

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

PROJECT="${PROJECT_ROOT}"
REF="${REF_FASTA}"

CHR2_BAM="ENA|OZ004312|OZ004312.1"
CHR2_VCF="ENA_OZ004312_OZ004312.1"
REGION="${CHR2_BAM}:11000000-16000000"

OUTDIR="${PROJECT}/cpl2a_analysis"
mkdir -p "${OUTDIR}/indels" "${OUTDIR}/delly" "${OUTDIR}/discordant" logs

source "${CONDA_SH}"
conda activate "${ENV_CPL2A}"
module load samtools 2>/dev/null || true
module load bcftools 2>/dev/null || true
command -v samtools >/dev/null || { echo "[ERROR] samtools not found" >&2; exit 1; }
command -v bcftools >/dev/null || { echo "[ERROR] bcftools not found" >&2; exit 1; }

BAM_ME1="${PROJECT}/samples/ME1/samfile/pip_1021/ME1.filtered.bam"
BAM_FL1="${PROJECT}/samples/FL1/samfile/pip_1021/FL1.filtered.bam"
BAM_CAL1="${PROJECT}/samples/CAL1/samfile/pip_1021/CAL1.filtered.bam"

for bam in "$BAM_ME1" "$BAM_FL1" "$BAM_CAL1"; do
    [[ -f "$bam" ]] || { echo "[ERROR] BAM not found: $bam" >&2; exit 1; }
done

# =========================================================================
# PART A: Indel calling with bcftools
# =========================================================================
echo "=== Part A: Indel calling ==="

INDEL_VCF="${OUTDIR}/indels/cpl2a_indels_all.vcf.gz"

# mpileup the three reference BAMs at the wide region, then call variants
# including indels. Pool-seq BAMs have high depth, so raise max-depth.
echo "  Running bcftools mpileup + call..."
bcftools mpileup \
    -q 10 \
    -Q 20 \
    -d 500 \
    -a FORMAT/AD,FORMAT/DP \
    -f "$REF" \
    -r "$REGION" \
    "$BAM_ME1" "$BAM_FL1" "$BAM_CAL1" \
    | bcftools call \
        -mv \
        --ploidy 1 \
        -Oz \
        -o "$INDEL_VCF"

bcftools index -f "$INDEL_VCF"

# Count total variants and indels
n_total=$(bcftools view -H "$INDEL_VCF" | wc -l)
n_indel=$(bcftools view -H -v indels "$INDEL_VCF" | wc -l)
n_snp=$(bcftools view -H -v snps "$INDEL_VCF" | wc -l)
echo "  Wide region: ${n_total} variants (${n_snp} SNPs, ${n_indel} indels)"

# Extract indels only in the candidate interval
echo "  Extracting candidate-interval indels..."
INDEL_CAND="${OUTDIR}/indels/cpl2a_indels_candidate.tsv"
{
    echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tTYPE\tME1_GT\tME1_DP\tFL1_GT\tFL1_DP\tCAL1_GT\tCAL1_DP"
    bcftools view -v indels -r "${CHR2_BAM}:13500000-14000000" "$INDEL_VCF" \
        | bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%TYPE[\t%GT\t%DP]\n'
} > "$INDEL_CAND"

n_cand=$(tail -n +2 "$INDEL_CAND" | wc -l)
echo "  Candidate interval: ${n_cand} indels"

# Find indels where PIP (ME1) and MOL (CAL1) are fixed-different
# GT field: 0 = ref, 1 = alt (haploid calling)
echo "  Finding PIP-MOL fixed-different indels..."
FIXED_DIFF="${OUTDIR}/indels/cpl2a_indels_pip_mol_fixed.tsv"
{
    head -1 "$INDEL_CAND"
    tail -n +2 "$INDEL_CAND" | awk -F'\t' '
        # ME1_GT=$7, CAL1_GT=$11
        # Fixed-different: one is 0, other is 1 (or vice versa)
        ($7 == "0" && $11 == "1") || ($7 == "1" && $11 == "0")
    '
} > "$FIXED_DIFF"

n_fixed=$(tail -n +2 "$FIXED_DIFF" | wc -l)
echo "  PIP-MOL fixed-different indels: ${n_fixed}"

if [[ "$n_fixed" -gt 0 ]]; then
    echo ""
    echo "  === Fixed PIP-MOL indel differences at CpL2a ==="
    head -30 "$FIXED_DIFF" | column -t -s $'\t' || true
fi

# =========================================================================
# PART B: Extract chr2 calls from existing DELLY output
# =========================================================================
echo ""
echo "=== Part B: DELLY chr2 extraction ==="

DELLY_DIR="${PROJECT}/discordant_pair_analysis/delly"

for SAMPLE in ME1 FL1 CAL1; do
    BCF="${DELLY_DIR}/${SAMPLE}_inversions.bcf"
    if [[ ! -f "$BCF" ]]; then
        echo "  [WARN] DELLY BCF not found for ${SAMPLE}: $BCF — skipping"
        continue
    fi

    # Extract all chr2 calls (any filter status)
    CHR2_OUT="${OUTDIR}/delly/${SAMPLE}_inversions_chr2_all.tsv"
    {
        echo -e "CHROM\tSTART\tEND\tSVLEN\tPE\tSR\tMAPQ\tCT\tFILTER"
        bcftools query \
            -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%INFO/PE\t%INFO/SR\t%INFO/MAPQ\t%INFO/CT\t%FILTER\n' \
            -r "${CHR2_BAM}" \
            "$BCF" 2>/dev/null || \
        bcftools query \
            -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%INFO/PE\t%INFO/SR\t%INFO/MAPQ\t%INFO/CT\t%FILTER\n' \
            "$BCF" | awk -v chr="$CHR2_BAM" '$1 == chr'
    } > "$CHR2_OUT"

    n_chr2=$(($(wc -l < "$CHR2_OUT") - 1))
    echo "  ${SAMPLE}: ${n_chr2} chr2 inversion calls"

    # Filter to CpL2a region (11-16 Mb, generous)
    CpL2a_OUT="${OUTDIR}/delly/${SAMPLE}_inversions_cpl2a.tsv"
    awk -F'\t' 'NR==1 || ($2 >= 11000000 && $2 <= 16000000)' "$CHR2_OUT" > "$CpL2a_OUT"
    n_cpl2a=$(($(wc -l < "$CpL2a_OUT") - 1))
    echo "  ${SAMPLE}: ${n_cpl2a} in CpL2a wide region"

    if [[ "$n_cpl2a" -gt 0 ]]; then
        head -20 "$CpL2a_OUT" | column -t -s $'\t' || true
    fi
    echo ""
done

# =========================================================================
# PART C: Discordant and split reads at CpL2a
# =========================================================================
echo "=== Part C: Discordant reads ==="

for SAMPLE in ME1 FL1 CAL1; do
    BAM="${PROJECT}/samples/${SAMPLE}/samfile/pip_1021/${SAMPLE}.filtered.bam"
    DISC="${OUTDIR}/discordant/${SAMPLE}_cpl2a_discordant.bam"

    echo "  ${SAMPLE}: extracting discordant pairs at CpL2a candidate..."
    # -F 1294: exclude unmapped, secondary, supplementary, duplicate, failed QC
    # -f 1: must be paired
    # Discordant pairs (mate on different chr, or insert size > 10kb) are
    # counted below
    samtools view -b -F 1294 -f 1 "$BAM" "${CHR2_BAM}:13500000-14000000" \
        | samtools sort -o "$DISC"
    samtools index "$DISC"

    n_disc=$(samtools view -c "$DISC")
    echo "  ${SAMPLE}: ${n_disc} reads in candidate interval (paired, primary)"

    # Count reads with mate on different chromosome
    n_inter=$(samtools view "$DISC" | awk '$7 != "=" && $7 != "*"' | wc -l)
    echo "  ${SAMPLE}: ${n_inter} with mate on different chromosome"

    # Count reads with insert size > 10kb (potential large indel/SV)
    n_large=$(samtools view "$DISC" | awk 'function abs(x) {return x<0?-x:x} abs($9) > 10000' | wc -l)
    echo "  ${SAMPLE}: ${n_large} with insert size > 10 kb"
    echo ""
done

echo "[✓] All CpL2a indel/SV analyses complete. Outputs in: $OUTDIR"
