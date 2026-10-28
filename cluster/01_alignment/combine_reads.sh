#!/bin/bash
# =============================================================================
# cluster/01_alignment/combine_reads.sh
# =============================================================================
# Concatenates each pool's adapter- and quality-trimmed reads from the two
# sequencing lanes (L001 then L004) into one R1 and one R2 file. Trimming
# (Trimmomatic v0.33) was done by the sequencing facility before delivery.
# Run from the project root: bash combine_reads.sh
#
# Step: 1, alignment (cluster/README.md)
# Manuscript: Supplementary Methods "Read processing, alignment, variant
#             calling and software versions"
#
# Inputs:  locations.txt (one pool per line)
#          <pool>/ILLUMINA_DATA/*L00{1,4}_R{1,2}_trimmed.fastq.gz (facility delivery)
# Outputs: <pool>/full_reads/full_R1.fastq.gz, full_R2.fastq.gz
# =============================================================================

# Read each location name from locations.txt
while read loc; do
    # Skip empty lines or commented lines
    [[ -z "$loc" || "$loc" == \#* ]] && continue

    illumina_dir="${loc}/ILLUMINA_DATA"
    outdir="${loc}/full_reads"
    mkdir -p "$outdir"

    # Find and sort matching files (L001 first, then L004)
    R1_files=$(echo "$illumina_dir"/*L001_R1_trimmed.fastq.gz "$illumina_dir"/*L004_R1_trimmed.fastq.gz 2>/dev/null)
    R2_files=$(echo "$illumina_dir"/*L001_R2_trimmed.fastq.gz "$illumina_dir"/*L004_R2_trimmed.fastq.gz 2>/dev/null)

    if [[ -n "$R1_files" && -n "$R2_files" ]]; then
        echo "Combining reads for $loc..."
        cat $R1_files > "$outdir/full_R1.fastq.gz"
        cat $R2_files > "$outdir/full_R2.fastq.gz"
    else
        echo "Skipping $loc: missing one or both R1/R2 read files"
    fi
done < locations.txt

