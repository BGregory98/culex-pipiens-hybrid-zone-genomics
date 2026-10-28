#!/bin/bash
# =============================================================================
# cluster/01_alignment/align_convert_sort.sh
# =============================================================================
# Aligns each pool's reads to idCulPipi1 with BWA-MEM2 2.2.1 (mem -k 19 -M,
# read group = pool name), then converts to a coordinate-sorted BAM, indexes
# it and writes flagstat with samtools 1.21. One array task per line of
# locations.txt (33 pools sequenced; 27 used downstream).
#
# Step: 1, alignment (cluster/README.md)
# Manuscript: Supplementary Methods "Read processing, alignment, variant
#             calling and software versions"; "Confirmation of inversion
#             breakpoints" (-k 19 -M)
#
# Inputs:  locations.txt
#          <pool>/full_reads/full_R{1,2}.fastq.gz (combine_reads.sh)
#          reference: REF_FASTA (BWA-MEM2 index files alongside)
# Outputs: <pool>/samfile/pip_1021/<pool>.bam (+ .bai), <pool>.flagstat.txt
# =============================================================================
#SBATCH --job-name=BWA_MEM2
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=16:00:00
#SBATCH -o ./sbatch_logs-BWAMEM2/%A_%a.out
#SBATCH -e ./sbatch_logs-BWAMEM2/%A_%a.err
#SBATCH --array=1-33

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

# --- Configuration ---
REF_INDEX="${REF_FASTA}"
SAMTOOLS="${SAMTOOLS_1_21_DIR}"

# Sample
LOCATION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" locations.txt)
R1="${LOCATION}/full_reads/full_R1.fastq.gz"
R2="${LOCATION}/full_reads/full_R2.fastq.gz"

# Output paths
SAMDIR="${LOCATION}/samfile/pip_1021"
SAMFILE="${SAMDIR}/${LOCATION}.sam"
BAMFILE="${SAMDIR}/${LOCATION}.bam"
BAMINDEX="${BAMFILE}.bai"
FLAGSTAT="${SAMDIR}/${LOCATION}.flagstat.txt"

# Create output directory
mkdir -p "$SAMDIR"

# Run alignment
if [[ -f "$R1" && -f "$R2" ]]; then
    echo "[$LOCATION] Running BWA-MEM2 alignment..."

    RG="@RG\\tID:${LOCATION}\\tSM:${LOCATION}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1"
    "$BWA_MEM2" mem -k 19 -t 8 -M -R "$RG" "$REF_INDEX" "$R1" "$R2" > "$SAMFILE"

    echo "[$LOCATION] Alignment complete. Converting to sorted BAM..."

    # Convert SAM to sorted BAM
    "$SAMTOOLS/samtools" view -bS "$SAMFILE" | "$SAMTOOLS/samtools" sort -o "$BAMFILE"

    if [ $? -eq 0 ]; then
        echo "[$LOCATION] BAM sorting complete. Cleaning up SAM..."
        rm "$SAMFILE"

        echo "[$LOCATION] Indexing BAM..."
        "$SAMTOOLS/samtools" index "$BAMFILE"

        echo "[$LOCATION] Generating flagstat report..."
        "$SAMTOOLS/samtools" flagstat "$BAMFILE" > "$FLAGSTAT"

        echo "[$LOCATION] All steps completed successfully."
    else
        echo "[$LOCATION] BAM sorting failed."
        exit 1
    fi
else
    echo "[$LOCATION] Missing FASTQ files. Skipping..."
    exit 1
fi

