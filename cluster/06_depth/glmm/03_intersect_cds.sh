#!/bin/bash
# cluster/06_depth/glmm/03_intersect_cds.sh
# =========================================
#
# Step 3 of the callable-depth pipeline (cluster/README.md step 8; Methods
# "Callable-position read-depth analysis"). Submitted through
# 03_intersect_cds.slurm.
#
# Input:  the callable_depth TSV from step 2
#         (chrom, pos, then 27 per-pool depths)
# Output: the same data filtered to CDS positions only, with a gene_id column
#         added between pos and the depth columns.
#
# The CDS annotation comes from the lifted Cx. quinquefasciatus GFF
# (reference_genomes/quinq_on_pip.gff3_polished). We take the union of all
# CDS features per gene_id, so a position falling in any CDS of any
# transcript of a gene is tagged with that gene's ID. Overlapping CDS ranges
# are merged with bedtools merge so a position is counted once per
# (pool, position) row. The merge does not group by gene: where CDS of two
# genes overlap or abut, the merged interval is labelled with both IDs
# (comma-separated) and counted as one unit downstream.
#
# Usage
# -----
#     bash 03_intersect_cds.sh \
#         <callable_depth.tsv.gz> \
#         <quinq_on_pip.gff3_polished> \
#         <output_dir>
#
# Outputs (all in <output_dir>):
#   cds_per_gene_merged.bed       BED4 (chrom, start, end, gene_id) of merged
#                                 CDS coordinates per gene_id. Reusable across
#                                 reruns; cache to avoid re-deriving.
#   callable_cds_depth.tsv.gz     The callable depth matrix restricted to CDS
#                                 positions, with gene_id column injected.
#
# Performance
# -----------
# The bedtools intersect step streams the callable file (which is sorted by
# chrom/pos already since step 2 writes in mpileup order). Wall-clock is
# minutes, not hours.

set -Eeuo pipefail
trap 'echo "[ERROR] line $LINENO status $?"' ERR

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <callable_depth.tsv.gz> <quinq_on_pip.gff3_polished> <output_dir>" >&2
    exit 1
fi

CALLABLE_DEPTH="$1"
GFF="$2"
OUTDIR="$3"

mkdir -p "$OUTDIR"

# Locate bedtools. Try module first, then conda env, then PATH.
if ! command -v bedtools >/dev/null 2>&1; then
    module load bedtools 2>/dev/null || true
fi
if ! command -v bedtools >/dev/null 2>&1; then
    echo "[ERROR] bedtools not in PATH and 'module load bedtools' failed" >&2
    exit 1
fi
echo "[INFO] using $(which bedtools) ($(bedtools --version))"

# ---------------------------------------------------------------------------
# Part A: build the merged-per-gene CDS BED from the GFF
# ---------------------------------------------------------------------------
#
# The GFF is GFF3 with CDS records having a Parent= attribute pointing at
# the mRNA, and the mRNA having a Parent= attribute pointing at the gene.
# Liftoff sometimes adds suffixes (_1, _2, ...) to gene IDs for multi-copy
# lifts; the awk below keeps them, so each lifted copy is its own gene_id.

CDS_BED="${OUTDIR}/cds_per_gene_merged.bed"

if [[ -f "$CDS_BED" ]]; then
    echo "[INFO] reusing cached CDS BED: $CDS_BED"
else
    echo "[INFO] building merged-per-gene CDS BED from $GFF"

    # The GFF may have CDS Parent= as transcript IDs that don't follow the
    # gene_id naming. To map CDS -> gene_id reliably, we extract via the
    # gene= attribute that liftoff adds, falling back to the locus_tag
    # convention used in the Cx. quinq annotation. The exact attribute key
    # may need adjusting after a head-check of the GFF. The awk below
    # prefers an explicit gene_id= attribute, then falls back to Parent=
    # with transcript-suffix stripping.

    # Step A1: pull gene-level coordinate records.
    # Step A2: filter to CDS, extract gene_id from the attributes column.
    # Step A3: convert to BED (0-based, end-exclusive), sort, merge by gene_id.

    awk -F '\t' '
        BEGIN { OFS = "\t" }
        $3 == "CDS" {
            # Attributes column is field 9. Parse semicolon-separated key=value pairs.
            # Preferred attribute keys for gene identity, in priority order:
            #   gene_id, gene, Parent (with transcript-suffix stripping)
            n = split($9, attrs, ";")
            gene_id = ""
            parent = ""
            gene_name = ""
            for (i = 1; i <= n; i++) {
                if (match(attrs[i], /^gene_id=/)) {
                    gene_id = substr(attrs[i], 9)
                } else if (match(attrs[i], /^gene=/)) {
                    gene_name = substr(attrs[i], 6)
                } else if (match(attrs[i], /^Parent=/)) {
                    parent = substr(attrs[i], 8)
                }
            }
            # Fallback chain: gene_id > gene > stripped Parent
            if (gene_id == "" && gene_name != "") gene_id = gene_name
            if (gene_id == "" && parent != "") {
                # Parent is typically a transcript ID; strip transcript-level
                # suffix to get gene ID. Lifted-quinq IDs look like
                # CQUJHB014704 (gene) and CQUJHB014704-RA (transcript), so
                # strip a trailing dash+letters.
                sub(/-R[A-Z]$/, "", parent)
                # Liftoff multi-copy suffix (_1, _2) is at the gene level;
                # do not strip here (we want CQUJHB014704_2 to remain
                # distinct from CQUJHB014704 if user wants per-copy
                # resolution).
                gene_id = parent
            }
            if (gene_id == "") next  # skip records we cannot attribute

            # GFF is 1-based, inclusive; BED is 0-based, end-exclusive.
            start_bed = $4 - 1
            end_bed   = $5
            print $1, start_bed, end_bed, gene_id
        }
    ' "$GFF" \
        | sort -k1,1 -k4,4 -k2,2n \
        > "${OUTDIR}/cds_all.bed4.unsorted_for_merge"

    # Merge overlapping/adjacent CDS intervals. bedtools merge requires
    # sorting by chrom then start; -c 4 -o distinct keeps the gene_id label,
    # or the comma-separated IDs when intervals of different genes merge.
    sort -k1,1 -k2,2n "${OUTDIR}/cds_all.bed4.unsorted_for_merge" \
        | bedtools merge -i - -c 4 -o distinct \
        > "$CDS_BED"

    rm -f "${OUTDIR}/cds_all.bed4.unsorted_for_merge"

    n_intervals=$(wc -l < "$CDS_BED")
    n_genes=$(cut -f4 "$CDS_BED" | sort -u | wc -l)
    echo "[INFO] built $CDS_BED: $n_intervals merged CDS intervals across $n_genes gene_ids"
fi

# ---------------------------------------------------------------------------
# Part B: intersect callable positions with the per-gene CDS BED
# ---------------------------------------------------------------------------
#
# The callable depth file is (chrom, pos, depth1..depth27). bedtools wants
# BED-style intervals on stdin; we convert pos to a 1-bp interval
# (pos-1, pos) on the fly. After intersection, bedtools writes the input
# columns followed by the BED4 columns from the CDS file; we then strip
# the redundant chrom/start/end and keep gene_id in column 3 of the output.

CALLABLE_CDS_OUT="${OUTDIR}/callable_cds_depth.tsv.gz"

echo "[INFO] intersecting callable positions with CDS BED"

# Header for the output: chrom, pos, gene_id, then 27 pool columns (we
# preserve column order from the input).
input_header=$(zcat "$CALLABLE_DEPTH" | { IFS= read -r line; cat > /dev/null; printf '%s\n' "$line"; })
# input_header is "chrom\tpos\t<pool1>\t<pool2>\t...\t<pool27>"
# Insert "gene_id" between pos and the pool columns.
out_header=$(awk -F '\t' 'BEGIN{OFS="\t"} {
    printf "%s\t%s\tgene_id", $1, $2
    for (i = 3; i <= NF; i++) printf "\t%s", $i
    printf "\n"
}' <<< "$input_header")

{
    echo "$out_header"

    # Skip the header line of the callable file, convert remaining rows to a
    # bedtools-friendly format (chrom, pos-1, pos, depths...), intersect with
    # the CDS BED, then reformat the output.
    zcat "$CALLABLE_DEPTH" \
        | tail -n +2 \
        | awk -F '\t' 'BEGIN{OFS="\t"} {
            start = $2 - 1
            printf "%s\t%d\t%d", $1, start, $2
            for (i = 3; i <= NF; i++) printf "\t%s", $i
            printf "\n"
        }' \
        | bedtools intersect \
            -a - \
            -b "$CDS_BED" \
            -wa -wb \
            -sorted \
        | awk -F '\t' '
            BEGIN { OFS = "\t" }
            {
                # Output columns from -wa -wb:
                #   1..N_a: full record from -a (chrom, start, end, depths)
                #   then N_a+1..N_a+N_b: full record from -b (chrom, start, end, gene_id)
                # We discard the duplicate chrom/start/end on the -b side and
                # the start on the -a side. We keep -a chrom, -a end (pos),
                # gene_id from -b, and the depth columns.
                #
                # N_a = 3 fixed columns (chrom, start, end) + n_pools depth columns
                # We do not know n_pools at awk time, but we know -b adds exactly
                # 4 columns (chrom, start, end, gene_id), so gene_id is $NF.
                na = NF - 4
                # Print chrom, pos (which is the -a end), gene_id, then depths.
                printf "%s\t%s\t%s", $1, $3, $NF
                for (i = 4; i <= na; i++) printf "\t%s", $i
                printf "\n"
            }
        '
} | gzip > "$CALLABLE_CDS_OUT"

n_callable_cds=$(zcat "$CALLABLE_CDS_OUT" | tail -n +2 | wc -l)
echo "[INFO] callable CDS positions: $n_callable_cds"
echo "[INFO] output: $CALLABLE_CDS_OUT"
echo "[OK] done."
