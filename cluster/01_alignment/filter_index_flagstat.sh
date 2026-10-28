#!/bin/bash
# =============================================================================
# cluster/01_alignment/filter_index_flagstat.sh
# =============================================================================
# Keeps properly paired reads with MAPQ >= 10 (samtools view -bh -f 2 -q 10),
# indexes the filtered BAM and writes flagstat. One array task per line of
# locations.txt. The filtered BAMs are the input to every pileup.
#
# Step: 2, duplicate removal and filtering (cluster/README.md)
# Manuscript: Supplementary Methods "Read processing, alignment, variant
#             calling and software versions"
#
# Inputs:  locations.txt
#          <pool>/samfile/pip_1021/<pool>.dedup.bam (remove_duplicates.sh)
# Outputs: <pool>/samfile/pip_1021/<pool>.filtered.bam (+ .bai),
#          <pool>.filtered.flagstat.txt
# =============================================================================
#SBATCH --job-name=FilterBAM
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH -o ./sbatch_logs-FilterBAM/%A_%a.out
#SBATCH -e ./sbatch_logs-FilterBAM/%A_%a.err
#SBATCH --array=1-33

# Load Samtools
module load samtools

# Identify sample
LOCATION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" locations.txt)

INFILE="${LOCATION}/samfile/pip_1021/${LOCATION}.dedup.bam"
OUTFILE="${LOCATION}/samfile/pip_1021/${LOCATION}.filtered.bam"
INDEXFILE="${OUTFILE}.bai"
FLAGSTAT="${LOCATION}/samfile/pip_1021/${LOCATION}.filtered.flagstat.txt"

echo "[$LOCATION] Starting filtering for proper pairs with MAPQ ≥10..."

if [ -f "$INFILE" ]; then
    # Filter BAM
    samtools view -bh -f 2 -q 10 "$INFILE" > "$OUTFILE"

    if [ $? -ne 0 ]; then
        echo "[$LOCATION] Filtering failed. Exiting..."
        exit 1
    fi

    echo "[$LOCATION] Filtering complete. Indexing..."

    # Index filtered BAM
    samtools index "$OUTFILE"

    if [ $? -ne 0 ]; then
        echo "[$LOCATION] Indexing failed."
    else
        echo "[$LOCATION] Index created: $INDEXFILE"
    fi

    # Compute flagstat
    echo "[$LOCATION] Computing flagstat..."
    samtools flagstat "$OUTFILE" > "$FLAGSTAT"
    echo "[$LOCATION] Flagstat written to $FLAGSTAT"

else
    echo "[$LOCATION] Input BAM not found. Skipping."
fi

