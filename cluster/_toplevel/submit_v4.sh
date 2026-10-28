#!/bin/bash
# =============================================================================
# cluster/_toplevel/submit_v4.sh
# =============================================================================
# Submits the pileup and variant-calling chain with afterok dependencies:
# mpileup_v4_array.sh (27 tasks) -> mpileup_v4_merge.sh -> PoolSNP_v4.sh.
# Run from the directory holding the three scripts (the project root in the
# published run): bash submit_v4.sh
#
# Step: 3-4, pileup and variant calling (cluster/README.md)
# Inputs:  see the three submitted scripts
# Outputs: SLURM job IDs printed to the terminal
# =============================================================================

set -Eeuo pipefail
cd "$(dirname "$(readlink -f "$0")")"

array_id=$(sbatch --parsable mpileup_v4_array.sh)
merge_id=$(sbatch --parsable --dependency=afterok:${array_id} mpileup_v4_merge.sh)
psp_id=$(sbatch   --parsable --dependency=afterok:${merge_id} PoolSNP_v4.sh)

echo "submitted:"
echo "  mpileup_v4_array  -> jobid ${array_id}  (27 parallel tasks)"
echo "  mpileup_v4_merge  -> jobid ${merge_id}  (depends on afterok:${array_id})"
echo "  PoolSNP_v4        -> jobid ${psp_id}    (depends on afterok:${merge_id})"
echo ""
echo "monitor:  squeue -u \$USER"
