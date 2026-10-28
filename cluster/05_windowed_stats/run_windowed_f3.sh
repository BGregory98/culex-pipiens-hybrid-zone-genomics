#!/bin/bash
# =============================================================================
# cluster/05_windowed_stats/run_windowed_f3.sh
# =============================================================================
# SLURM array driver for windowed f3*(focal; ref1, ref2): array task i takes
# the i-th field pool (sorted names from data/metadata/sample_categories.csv)
# as the focal population and runs windowed_f3.R. One submission per
# reference pairing, from the project root:
#   sbatch --array=1-24 --export=REF1=FL1,REF2=CAL1,WIN_SIZE=1000000 scripts/05_windowed_stats/run_windowed_f3.sh
#   sbatch --array=1-24 --export=REF1=ME1,REF2=CAL1,WIN_SIZE=1000000 scripts/05_windowed_stats/run_windowed_f3.sh
#   sbatch --array=1-24 --export=REF1=ME1,REF2=FL1,WIN_SIZE=1000000  scripts/05_windowed_stats/run_windowed_f3.sh
# The first task writes the field-pool list to logs/.field_pools_cache.txt;
# a task that reads the cache while it is being written stops with "no focal
# pool" and was resubmitted on its own (--array=<task>).
#
# Step: 7, windowed f3 (cluster/README.md)
# Manuscript: Methods "Windowed f3 admixture scan"
#
# Inputs:  SequencingData/data/pooldata_v4.rds,
#          SequencingData/data/windows_tbl_v4_<size>.rds (prepare_pooldata.R)
# Outputs: SequencingData/data/windowed_f3_parts_v4/
#            windowed_F3_<focal>_<ref1>_<ref2>_<size>.csv
# =============================================================================
#SBATCH --job-name=winF3
#SBATCH --output=logs/winF3_%A_%a.out
#SBATCH --error=logs/winF3_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -Eeuo pipefail

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs SequencingData/data/windowed_f3_parts_v4

module purge
module load R/4.5.3

REF1=${REF1:?Set REF1 in --export, e.g. REF1=FL1}
REF2=${REF2:?Set REF2 in --export, e.g. REF2=CAL1}
WIN_SIZE=${WIN_SIZE:-1000000}

# Format size label for filename
if   (( WIN_SIZE % 1000000 == 0 )); then SIZE_LABEL="$((WIN_SIZE / 1000000))M"
elif (( WIN_SIZE % 1000    == 0 )); then SIZE_LABEL="$((WIN_SIZE / 1000))K"
else                                     SIZE_LABEL="$WIN_SIZE"
fi

# Sorted field-pool list from sample_categories.csv, the same list the BayPass
# scripts use.
FIELD_POOLS_FILE="logs/.field_pools_cache.txt"
if [[ ! -s "$FIELD_POOLS_FILE" ]]; then
  Rscript -e '
    si <- read.csv("data/metadata/sample_categories.csv")
    fp <- sort(si$sample[si$sample_type == "field"])
    writeLines(fp, "logs/.field_pools_cache.txt")
  '
fi

mapfile -t FIELD_POPS < "$FIELD_POOLS_FILE"
FOCAL_POP="${FIELD_POPS[$SLURM_ARRAY_TASK_ID-1]:-}"

if [[ -z "$FOCAL_POP" ]]; then
  echo "ERROR: no focal pool for array task ${SLURM_ARRAY_TASK_ID}" >&2
  echo "Field pool list:" >&2
  printf '  %s\n' "${FIELD_POPS[@]}" >&2
  exit 1
fi

OUT_CSV="SequencingData/data/windowed_f3_parts_v4/windowed_F3_${FOCAL_POP}_${REF1}_${REF2}_${SIZE_LABEL}.csv"

echo "=== windowed f3* v4 ==="
echo "  Focal:     $FOCAL_POP"
echo "  Ref1:      $REF1"
echo "  Ref2:      $REF2"
echo "  Win size:  $WIN_SIZE bp ($SIZE_LABEL)"
echo "  Output:    $OUT_CSV"
echo "  Job:       $SLURM_JOB_ID  Array task: $SLURM_ARRAY_TASK_ID"
echo "  Host:      $(hostname)"
echo "  Started:   $(date)"

Rscript scripts/05_windowed_stats/windowed_f3.R "$FOCAL_POP" "$REF1" "$REF2" "$WIN_SIZE" "$OUT_CSV"

echo "  Finished:  $(date)"
