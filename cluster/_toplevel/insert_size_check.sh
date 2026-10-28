#!/bin/bash
# =============================================================================
# cluster/_toplevel/insert_size_check.sh
# =============================================================================
# Insert-size summary (n, mean, median, 0.1/1/99/99.9th percentiles) for the
# three reference pools from a 0.1% read subsample (samtools view -s 0.001) of
# the duplicate-removed BAMs, counting pairs with 0 < TLEN < 2,000 bp. Used to
# check that one set of DELLY parameters suits all three pools.
#
# Step: 10, inversion breakpoints (cluster/README.md)
# Manuscript: Supplementary Methods "Confirmation of inversion breakpoints"
#             (mean 220-246 bp; 99th percentile 712-781 bp)
#
# Inputs:  {CAL1,ME1,FL1}/samfile/pip_1021/<pool>.dedup.bam (01_alignment/remove_duplicates.sh)
# Outputs: summary lines in sbatch_logs-InsertSize/<jobid>.out
# =============================================================================
#SBATCH --job-name=InsertSize
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH -o ./sbatch_logs-InsertSize/%j.out
#SBATCH -e ./sbatch_logs-InsertSize/%j.err

source "${CLUSTER_CONFIG:-scripts/config.sh}" || exit 1

module load samtools

PROJECT="${PROJECT_ROOT}"

CAL1_BAM="${PROJECT}/CAL1/samfile/pip_1021/CAL1.dedup.bam"
ME1_BAM="${PROJECT}/ME1/samfile/pip_1021/ME1.dedup.bam"
FL1_BAM="${PROJECT}/FL1/samfile/pip_1021/FL1.dedup.bam"

mkdir -p ./sbatch_logs-InsertSize

for SAMPLE in CAL1 ME1 FL1; do
    case $SAMPLE in
        CAL1) BAM=$CAL1_BAM ;;
        ME1)  BAM=$ME1_BAM ;;
        FL1)  BAM=$FL1_BAM ;;
    esac

    echo "=== $SAMPLE ==="
    echo "BAM: $BAM"

    if [ ! -f "$BAM" ]; then
        echo "  MISSING: $BAM"
        continue
    fi

    samtools view -s 0.001 "$BAM" | \
        awk '$9>0 && $9<2000 {print $9}' | \
        sort -n | \
        awk 'BEGIN{c=0}
             {a[c++]=$1; sum+=$1}
             END{
               printf "  n=%d  mean=%.1f  median=%.1f  p1=%.1f  p99=%.1f  p0.1=%.1f  p99.9=%.1f\n",
                 c, sum/c, a[int(c/2)], a[int(c*0.01)], a[int(c*0.99)],
                 a[int(c*0.001)], a[int(c*0.999)]
             }'

    echo ""
done

echo "Done."
