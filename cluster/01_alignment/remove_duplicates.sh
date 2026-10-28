#!/bin/bash
# =============================================================================
# cluster/01_alignment/remove_duplicates.sh
# =============================================================================
# Removes PCR and optical duplicates with Picard 3.4.0 MarkDuplicates
# (REMOVE_DUPLICATES=true) and indexes the result. One array task per line of
# locations.txt. These BAMs, taken before the proper-pair filter, are also
# the DELLY input.
#
# Step: 2, duplicate removal and filtering (cluster/README.md)
# Manuscript: Supplementary Methods "Read processing, alignment, variant
#             calling and software versions"; "Confirmation of inversion
#             breakpoints"
#
# Inputs:  locations.txt
#          <pool>/samfile/pip_1021/<pool>.bam (align_convert_sort.sh)
# Outputs: <pool>/samfile/pip_1021/<pool>.dedup.bam (+ index),
#          <pool>.dedup.metrics.txt
# =============================================================================
#SBATCH --job-name=DedupBAM
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH -o ./sbatch_logs-DedupBAM/%A_%a.out
#SBATCH -e ./sbatch_logs-DedupBAM/%A_%a.err
#SBATCH --array=1-33

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

JAVA_HOME="${JAVA_HOME_DIR}"
JAVA_EXEC="$JAVA_HOME/bin/java"

# Identify sample
LOCATION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${PROJECT_ROOT}/locations.txt)
INPUT_BAM="${PROJECT_ROOT}/${LOCATION}/samfile/pip_1021/${LOCATION}.bam"
DEDUP_BAM="${PROJECT_ROOT}/${LOCATION}/samfile/pip_1021/${LOCATION}.dedup.bam"
METRICS="${PROJECT_ROOT}/${LOCATION}/samfile/pip_1021/${LOCATION}.dedup.metrics.txt"

echo "[$LOCATION] Starting MarkDuplicates"
echo "[$LOCATION] Java path: $JAVA_EXEC"
$JAVA_EXEC -version

# Check input BAM exists
if [ -f "$INPUT_BAM" ]; then
    echo "[$LOCATION] Input BAM found. Running Picard..."

    $JAVA_EXEC -Xmx12G -jar "$PICARD_JAR" MarkDuplicates \
        I="$INPUT_BAM" \
        O="$DEDUP_BAM" \
        M="$METRICS" \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT

    if [ $? -eq 0 ]; then
        echo "[$LOCATION] Duplicate removal complete."
    else
        echo "[$LOCATION] Picard encountered an error during duplicate removal."
    fi
else
    echo "[$LOCATION] Input BAM not found at $INPUT_BAM. Skipping."
fi

