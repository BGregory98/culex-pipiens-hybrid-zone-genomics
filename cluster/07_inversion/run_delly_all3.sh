#!/bin/bash
# =============================================================================
# cluster/07_inversion/run_delly_all3.sh
# =============================================================================
# DELLY 1.7.3 inversion calling (call -t INV -q 20) on the duplicate-removed,
# not proper-pair-filtered BAMs of the three reference pools (array task 1
# FL1/QUI, 2 CAL1/MOL, 3 ME1/PIP), then tiers of filtered calls:
#   - all calls
#   - PASS calls
#   - high-confidence PASS calls: PE >= 5, INFO/MAPQ >= 30, size 50 kb - 20 Mb
#     (removes cross-species micro-inversions, repeat artifacts and
#     chromosome-spanning false positives)
#   - PASS calls with both PE and SR support
#   - the chr3 window 148-163 Mb at the high-confidence and PE+SR tiers
# CAL1 and ME1 are within the Cx. pipiens complex, so cross-species mapping
# noise is lower than for FL1 (Cx. quinquefasciatus); FL1 runs through the
# same filters for consistency. The CAL1 high-confidence chr3 calls define
# the CpL3a breakpoints (149,402,438-161,021,387; R/helpers/loci.R).
#
# Step: 10, inversion breakpoints (cluster/README.md)
# Manuscript: Methods "Inversion detection and breakpoint confirmation";
#             Supplementary Methods "Confirmation of inversion breakpoints";
#             Supplementary Results "Inversion breakpoint coordinates"
#
# Inputs:  {FL1,CAL1,ME1}/samfile/pip_1021/<pool>.dedup.bam (01_alignment/remove_duplicates.sh)
#          reference: REF_FASTA
# Outputs: discordant_pair_analysis/delly/<pool>_inversions.bcf (+ .csi),
#          <pool>_inversions.{all,pass,highconf,pe_sr}.tsv,
#          <pool>_inversions.{highconf,pe_sr}.chr3_candidate.tsv
# =============================================================================
#SBATCH --job-name=DellyAll
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH -o ./sbatch_logs-Delly/%A_%a.out
#SBATCH -e ./sbatch_logs-Delly/%A_%a.err
#SBATCH --array=1-3

set -euo pipefail

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

module load samtools
module load bcftools 2>/dev/null || true

# DELLY binary: DELLY in cluster/config.sh unless DELLY_BIN is set.
DELLY_BIN="${DELLY_BIN:-${DELLY}}"

if [[ ! -x "$DELLY_BIN" ]]; then
    echo "ERROR: DELLY binary not found or not executable: $DELLY_BIN" >&2
    exit 1
fi

PROJECT="${PROJECT_ROOT}"
REF="${REF_FASTA}"
OUTDIR="${PROJECT}/discordant_pair_analysis/delly"
mkdir -p "$OUTDIR" ./sbatch_logs-Delly

# Pick sample by array index
case "$SLURM_ARRAY_TASK_ID" in
    1) SAMPLE="FL1"  ; BAM="${PROJECT}/FL1/samfile/pip_1021/FL1.dedup.bam" ;;
    2) SAMPLE="CAL1" ; BAM="${PROJECT}/CAL1/samfile/pip_1021/CAL1.dedup.bam" ;;
    3) SAMPLE="ME1"  ; BAM="${PROJECT}/ME1/samfile/pip_1021/ME1.dedup.bam" ;;
    *) echo "ERROR: invalid SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2 ; exit 1 ;;
esac

# Sanity-check inputs
[[ -f "$BAM" ]]      || { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }
[[ -f "${REF}.fai" ]] || { echo "[setup] Indexing reference..."; samtools faidx "$REF"; }

echo "=== DELLY inversion call: $SAMPLE ==="
echo "BAM:       $BAM"
echo "Reference: $REF"
echo "Started:   $(date)"

BCF="$OUTDIR/${SAMPLE}_inversions.bcf"

# ------- DELLY call (inversions only) -------
"$DELLY_BIN" call \
    -t INV \
    -g "$REF" \
    -o "$BCF" \
    -q 20 \
    "$BAM"

bcftools index -f "$BCF"
echo "Finished call: $(date)"

# ------- Tier 0: full call set (archival) -------
SUMMARY_ALL="$OUTDIR/${SAMPLE}_inversions.all.tsv"
{
    echo -e "CHROM\tSTART\tEND\tSVLEN\tPE\tSR\tMAPQ\tCT\tFILTER"
    bcftools query \
        -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%INFO/PE\t%INFO/SR\t%INFO/MAPQ\t%INFO/CT\t%FILTER\n' \
        "$BCF"
} > "$SUMMARY_ALL"

# ------- Tier 1: PASS only -------
SUMMARY_PASS="$OUTDIR/${SAMPLE}_inversions.pass.tsv"
awk -F'\t' 'NR==1 || $9 == "PASS"' "$SUMMARY_ALL" > "$SUMMARY_PASS"

# ------- Tier 2: high-confidence PASS -------
# PE >= 5, MAPQ >= 30, 50 kb <= size <= 20 Mb.
# Note: SR is "." when missing, so the support filter requires PE numeric only.
SUMMARY_HIGHCONF="$OUTDIR/${SAMPLE}_inversions.highconf.tsv"
awk -F'\t' 'NR==1 {print; next}
            $9 == "PASS" && $5 != "." && $5+0 >= 5 && $7+0 >= 30 {
                size = $3 - $2
                if (size >= 50000 && size <= 20000000) print
            }' "$SUMMARY_ALL" > "$SUMMARY_HIGHCONF"

# ------- Tier 3: PE + SR PASS (most reliable) -------
SUMMARY_PESR="$OUTDIR/${SAMPLE}_inversions.pe_sr.tsv"
awk -F'\t' 'NR==1 || ($9 == "PASS" && $5 != "." && $6 != ".")' \
    "$SUMMARY_ALL" > "$SUMMARY_PESR"

# ------- chr3 candidate region (148-163 Mb) at each tier -------
CHR3_REGEX='OZ004313'
for TIER in highconf pe_sr; do
    SRC="$OUTDIR/${SAMPLE}_inversions.${TIER}.tsv"
    DST="$OUTDIR/${SAMPLE}_inversions.${TIER}.chr3_candidate.tsv"
    awk -F'\t' -v re="$CHR3_REGEX" '
        NR==1 {print; next}
        $1 ~ re && $2 >= 148000000 && $2 <= 163000000 {print}
    ' "$SRC" > "$DST"
done

# ------- Report counts -------
echo ""
echo "=== Filter-tier counts for $SAMPLE ==="
printf "  All calls:                %d\n" $(($(wc -l < "$SUMMARY_ALL") - 1))
printf "  PASS:                     %d\n" $(($(wc -l < "$SUMMARY_PASS") - 1))
printf "  PASS + high-conf:         %d\n" $(($(wc -l < "$SUMMARY_HIGHCONF") - 1))
printf "  PASS + PE+SR:             %d\n" $(($(wc -l < "$SUMMARY_PESR") - 1))
printf "  chr3 148-163 Mb high-conf: %d\n" $(($(wc -l < "$OUTDIR/${SAMPLE}_inversions.highconf.chr3_candidate.tsv") - 1))
printf "  chr3 148-163 Mb PE+SR:     %d\n" $(($(wc -l < "$OUTDIR/${SAMPLE}_inversions.pe_sr.chr3_candidate.tsv") - 1))

echo ""
echo "=== chr3 candidate region (148-163 Mb), high-conf calls ==="
column -t -s $'\t' "$OUTDIR/${SAMPLE}_inversions.highconf.chr3_candidate.tsv" | head -30 || true

echo ""
echo "=== Outputs in $OUTDIR ==="
ls -la "$OUTDIR" | grep "^.*${SAMPLE}_inversions" || true
