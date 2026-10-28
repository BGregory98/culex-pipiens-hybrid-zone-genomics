#!/bin/bash
# =============================================================================
# cluster/06_depth/glmm/04b_aggregate_cds_depth_all27.slurm
# =============================================================================
# Reference-pool depth control for the CpL1a CYP6 cluster: runs step 4 again
# with --pool-set all, so the three reference pools (FL1/QUI, ME1/PIP,
# CAL1/MOL) are kept rather than dropped. Steps 1-3 already carry all 27
# pools, so they are not rerun. Writes *_all27.tsv beside the field-pool
# outputs without touching them.
#
# Step: 8, callable depth (cluster/README.md)
# Manuscript: Supplementary Methods "Callable-position read-depth analysis"
#             (reference-pool paragraph); Supplementary Results "Read depth in
#             the reference pools, and the haplotype that carries the extra
#             copies"
#
# Inputs:  data/glmm_v4/callable_cds_depth.tsv.gz, callable_depth_v4.tsv.gz
#          (steps 2-3), data/metadata/sample_categories.csv
# Outputs: data/glmm_v4/cds_depth_per_gene_pool_all27.tsv,
#          data/glmm_v4/genome_callable_depth_per_pool_all27.tsv (read by
#          R/08_loci/09_cyp6_reference_pool_depth.R)
#
# PATHS
# -----
# Each input path is taken from an environment variable if set, else from a
# list of candidate locations, else from a depth-limited search that refuses
# to choose between several copies of the same file name. The resolved paths
# are printed before any work starts.
#
# To pin any of them explicitly:
#
#   export CULEX_PROJ=$PROJECT_ROOT
#   export CULEX_GLMM_DIR=$CULEX_PROJ/scripts/06_depth/glmm
#   export CULEX_CALLABLE_CDS=$CULEX_PROJ/data/glmm_v4/callable_cds_depth.tsv.gz
#   export CULEX_CALLABLE_ALL=$CULEX_PROJ/data/glmm_v4/callable_depth_v4.tsv.gz
#   export CULEX_SAMPLE_CATEGORIES=$CULEX_PROJ/data/metadata/sample_categories.csv
#   export CULEX_OUT_DIR=$CULEX_PROJ/data/glmm_v4
#   sbatch 04b_aggregate_cds_depth_all27.slurm
#
# Or run it with --dry-run to resolve and report without doing any work:
#
#   bash 04b_aggregate_cds_depth_all27.slurm --dry-run
# =============================================================================
#SBATCH --job-name=glmm_04b_all27
#SBATCH --output=logs/glmm_04b_all27_%j.out
#SBATCH --error=logs/glmm_04b_all27_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=1

set -Eeuo pipefail

DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

PROJ="${CULEX_PROJ:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "$PROJ"
echo "[INFO] project root: $PROJ"

# --- resolution helpers ------------------------------------------------------
# find_file <name> <candidate dirs...> : take the first candidate directory that
# has the file. Only if none does, fall back to a depth-limited search -- and if
# that search finds MORE THAN ONE copy, refuse to choose.
#
# Refusing matters: an older file of the same name and column count would pass
# every header check below while producing wrong numbers, so an ambiguous match
# is left to the operator rather than to find's traversal order.
find_file() {
  local name="$1"; shift
  local d
  for d in "$@"; do
    if [[ -s "$d/$name" ]]; then echo "$d/$name"; return 0; fi
  done
  local -a hits=()
  while IFS= read -r line; do hits+=("${line#./}"); done < <(
    find . -maxdepth 5 -type f -name "$name" \
      -not -path './.git/*' -not -path './logs/*' 2>/dev/null | sort)
  case ${#hits[@]} in
    0) return 1 ;;
    1) echo "${hits[0]}"; return 0 ;;
    *) echo "[ERROR] $name is ambiguous -- ${#hits[@]} copies found and none in a" >&2
       echo "[ERROR] known location, so refusing to guess. Candidates:" >&2
       printf '[ERROR]     %s\n' "${hits[@]}" >&2
       echo "[ERROR] Pin the right one explicitly, e.g." >&2
       echo "[ERROR]     export CULEX_DATA_DIR=<dir holding the v4 step-2/3 outputs>" >&2
       echo "[ERROR] or set the per-file variable (see the header of this script)." >&2
       return 2 ;;
  esac
}

resolve() {  # resolve <label> <env value> <name> <candidate dirs...>
  local label="$1" envval="$2" name="$3"; shift 3
  local path
  if [[ -n "$envval" ]]; then
    path="$envval"
    if [[ ! -s "$path" ]]; then
      echo "[ERROR] $label pinned by environment but missing or empty: $path" >&2
      return 1
    fi
  else
    local rc=0
    path=$(find_file "$name" "$@") || rc=$?
    if [[ $rc -eq 2 ]]; then
      echo "[ERROR] (while resolving $label)" >&2
      return 1
    elif [[ $rc -ne 0 ]]; then
      echo "[ERROR] could not locate $label ($name)." >&2
      echo "[ERROR]   looked in: $*" >&2
      echo "[ERROR]   then searched $PROJ to depth 5 and found nothing." >&2
      return 1
    fi
  fi
  printf '%-22s %s\n' "[OK] $label" "$path"
  echo "$path" > "$RESOLVED_TMP"
}

RESOLVED_TMP=$(mktemp)
trap 'rm -f "$RESOLVED_TMP"' EXIT
get() { resolve "$@" >&2 && cat "$RESOLVED_TMP"; }

echo "[INFO] resolving inputs ..."

# The R script lives beside this one, wherever "beside" now is.
R_SCRIPT=$(get "step-4 R script" "${CULEX_R_SCRIPT:-}" "04_aggregate_cds_depth.R" \
  "${CULEX_GLMM_DIR:-/nonexistent}" \
  "scripts/06_depth/glmm" "scripts/05_depth/glmm" "cluster/06_depth/glmm" \
  "scripts/glmm" "scripts/06_depth" "scripts")

# Step-3 output: chrom, pos, gene_id, then one column per pool.
CALLABLE_CDS=$(get "callable CDS depth" "${CULEX_CALLABLE_CDS:-}" "callable_cds_depth.tsv.gz" \
  "${CULEX_DATA_DIR:-/nonexistent}" \
  "data/glmm_v4" "data/glmm" "glmm_v4" "data")

# Step-2 output: chrom, pos, then one column per pool. Needed for the
# genome-wide callable-depth denominator.
CALLABLE_ALL=$(get "callable all depth" "${CULEX_CALLABLE_ALL:-}" "callable_depth_v4.tsv.gz" \
  "${CULEX_DATA_DIR:-/nonexistent}" \
  "data/glmm_v4" "data/glmm" "glmm_v4" "data")

# Metadata: sample, sample_type, sample_type_2.
SAMPLE_CATS=$(get "sample categories" "${CULEX_SAMPLE_CATEGORIES:-}" "sample_categories.csv" \
  "${CULEX_META_DIR:-/nonexistent}" \
  "data/metadata" "data" "metadata")

# Outputs go beside the step-2/3 products unless told otherwise.
OUT_DIR="${CULEX_OUT_DIR:-$(dirname "$CALLABLE_CDS")}"
printf '%-22s %s\n' "[OK] output dir" "$OUT_DIR"

echo
echo "[INFO] this job will write:"
echo "[INFO]   $OUT_DIR/cds_depth_per_gene_pool_all27.tsv"
echo "[INFO]   $OUT_DIR/genome_callable_depth_per_pool_all27.tsv"
echo "[INFO] it will NOT touch the GLMM inputs cds_depth_per_gene_pool.tsv or"
echo "[INFO] genome_callable_depth_per_pool.tsv."
echo

# --- sanity-check the inputs actually look right ------------------------------
# Cheap header checks: the reference columns must be present before the
# 90-minute run starts.
echo "[INFO] checking headers ..."

# NB: `zcat file | head -1` exits 141 (SIGPIPE) under `set -o pipefail`, because
# head closes the pipe while the decompressor is still writing. With `set -e`
# that aborts the job. pipefail is disabled inside this subshell only, so the
# pipeline reports head's status instead.
read_header() {
  local f="$1"
  ( set +o pipefail
    if [[ "$f" == *.gz ]]; then gzip -dc -- "$f" | head -1; else head -1 -- "$f"; fi )
}
CDS_HDR=$(read_header "$CALLABLE_CDS")
ALL_HDR=$(read_header "$CALLABLE_ALL")
if [[ -z "$CDS_HDR" || -z "$ALL_HDR" ]]; then
  echo "[ERROR] could not read a header line from one of the inputs." >&2
  echo "[ERROR]   $CALLABLE_CDS" >&2
  echo "[ERROR]   $CALLABLE_ALL" >&2
  exit 1
fi
n_cds_cols=$(awk -F'\t' '{print NF; exit}' <<<"$CDS_HDR")
n_all_cols=$(awk -F'\t' '{print NF; exit}' <<<"$ALL_HDR")
echo "[INFO]   callable CDS: $n_cds_cols columns (expect 3 + n_pools)"
echo "[INFO]   callable all: $n_all_cols columns (expect 2 + n_pools)"

missing_ref=0
for ref in FL1 ME1 CAL1; do
  if grep -qw -- "$ref" <<<"$CDS_HDR" && grep -qw -- "$ref" <<<"$ALL_HDR"; then
    echo "[INFO]   reference pool $ref: present in both"
  else
    echo "[ERROR]  reference pool $ref: NOT present in both headers" >&2
    missing_ref=1
  fi
done
if [[ $missing_ref -ne 0 ]]; then
  echo "[ERROR] The reference pools are the entire point of this run. Their absence" >&2
  echo "[ERROR] means these files came from a field-only locations list, so steps 2" >&2
  echo "[ERROR] and 3 must be re-run against the 27-pool mpileup first (step 2 is" >&2
  echo "[ERROR] the expensive one: 3-day walltime, 16G)." >&2
  exit 1
fi

if [[ $DRY_RUN -eq 1 ]]; then
  echo
  echo "[OK] --dry-run: everything resolved and the headers look right. Nothing run."
  exit 0
fi

mkdir -p "$OUT_DIR" logs
module load R/4.5.3 || echo "[WARN] 'module load R/4.5.3' failed; relying on Rscript already on PATH"

Rscript "$R_SCRIPT" \
  --callable-cds      "$CALLABLE_CDS" \
  --callable-all      "$CALLABLE_ALL" \
  --sample-categories "$SAMPLE_CATS" \
  --min-positions     50 \
  --pool-set          all \
  --out-dir           "$OUT_DIR"

echo
echo "[OK] outputs for the R stage (data/glmm_v4/):"
echo "     $OUT_DIR/cds_depth_per_gene_pool_all27.tsv"
echo "     $OUT_DIR/genome_callable_depth_per_pool_all27.tsv"
echo "[OK] next: R/08_loci/09_cyp6_reference_pool_depth.R"
