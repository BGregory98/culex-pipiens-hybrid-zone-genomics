#!/usr/bin/env python3
# Extracts the gene, mRNA and CDS records of the target genes (the six CpL1a
# CYP6 paralogs and the P450 panel) from a GFF3, and writes a subset GFF3, the
# transcript and gene ID lists and a per-gene summary table.
# Step 12 (cluster/README.md). Run as:
#   extract_target_models.py --ids <gene_ids.txt> --gff <annotation.gff3> --out-prefix <prefix>
import argparse
from collections import defaultdict

def parse_attrs(s):
    out = {}
    for item in s.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out

def uniq_attr(attr_list, key):
    vals = []
    for d in attr_list:
        if key in d and d[key] != "":
            vals.append(d[key])
    return ",".join(sorted(set(vals)))

ap = argparse.ArgumentParser()
ap.add_argument("--ids", required=True, help="file with target gene IDs, one per line")
ap.add_argument("--gff", required=True, help="input GFF3/GFF")
ap.add_argument("--out-prefix", required=True, help="output prefix")
args = ap.parse_args()

targets = {x.strip() for x in open(args.ids) if x.strip()}
gene_like = {"protein_coding_gene", "gene"}
mrna_like = {"mRNA", "transcript"}

gene_rows = defaultdict(list)
mrna_rows = defaultdict(list)
gene_attrs = defaultdict(list)
mrna_attrs = defaultdict(list)
gene_coords = defaultdict(list)
gene_to_tids = defaultdict(list)
tid_to_gene = {}

# pass 1: gene + mRNA discovery
with open(args.gff) as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        seqid, source, feature, start, end, score, strand, phase, attr_s = cols
        attrs = parse_attrs(attr_s)

        if feature in gene_like:
            gid = attrs.get("ID")
            if gid in targets:
                gene_rows[gid].append(line)
                gene_attrs[gid].append(attrs)
                gene_coords[gid].append(f"{seqid}:{start}-{end}({strand})")

        elif feature in mrna_like:
            parents = [p for p in attrs.get("Parent", "").split(",") if p]
            hit_genes = [p for p in parents if p in targets]
            if hit_genes:
                tid = attrs.get("ID")
                if tid:
                    for gid in hit_genes:
                        gene_to_tids[gid].append(tid)
                        tid_to_gene[tid] = gid
                        mrna_rows[gid].append(line)
                        mrna_attrs[gid].append(attrs)

target_tids = set(tid_to_gene.keys())

# pass 2: write subset GFF in original order
subset_gff = args.out_prefix + ".subset.gff3"
with open(args.gff) as fh, open(subset_gff, "w") as out:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        feature = cols[2]
        attrs = parse_attrs(cols[8])
        parents = [p for p in attrs.get("Parent", "").split(",") if p]
        rec_id = attrs.get("ID", "")

        keep = False
        if feature in gene_like and rec_id in targets:
            keep = True
        elif feature in mrna_like and any(p in targets for p in parents):
            keep = True
        elif feature == "CDS" and any(p in target_tids for p in parents):
            keep = True

        if keep:
            out.write(line)

# transcript ID list
with open(args.out_prefix + ".transcript_ids.txt", "w") as out:
    for gid in sorted(targets):
        tids = sorted(set(gene_to_tids.get(gid, [])))
        for tid in tids:
            out.write(tid + "\n")

# gene ID list actually found
with open(args.out_prefix + ".gene_ids_found.txt", "w") as out:
    for gid in sorted(targets):
        if gid in gene_rows:
            out.write(gid + "\n")

# summary table
with open(args.out_prefix + ".summary.tsv", "w") as out:
    header = [
        "gene_id",
        "gene_feature_n",
        "mrna_feature_n",
        "transcript_ids",
        "coords",
        "description",
        "coverage",
        "sequence_ID",
        "valid_ORFs",
        "extra_copy_number",
        "copy_num_ID",
        "mrna_valid_ORF",
        "mrna_matches_ref_protein",
        "mrna_missing_stop_codon"
    ]
    out.write("\t".join(header) + "\n")
    for gid in sorted(targets):
        out.write("\t".join([
            gid,
            str(len(gene_rows.get(gid, []))),
            str(len(mrna_rows.get(gid, []))),
            ",".join(sorted(set(gene_to_tids.get(gid, [])))),
            ",".join(gene_coords.get(gid, [])),
            uniq_attr(gene_attrs.get(gid, []), "description"),
            uniq_attr(gene_attrs.get(gid, []), "coverage"),
            uniq_attr(gene_attrs.get(gid, []), "sequence_ID"),
            uniq_attr(gene_attrs.get(gid, []), "valid_ORFs"),
            uniq_attr(gene_attrs.get(gid, []), "extra_copy_number"),
            uniq_attr(gene_attrs.get(gid, []), "copy_num_ID"),
            uniq_attr(mrna_attrs.get(gid, []), "valid_ORF"),
            uniq_attr(mrna_attrs.get(gid, []), "matches_ref_protein"),
            uniq_attr(mrna_attrs.get(gid, []), "missing_stop_codon"),
        ]) + "\n")
