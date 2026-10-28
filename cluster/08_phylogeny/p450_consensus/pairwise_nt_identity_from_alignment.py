#!/usr/bin/env python3
# As pairwise_identity_from_alignment.py, for the whole-gene nucleotide
# alignments (<gene>.gene.aln.fa), skipping gaps and N. Saved as
# data/consensus/p450s/pairwise_gene_nt_identity.tsv (read by
# R/08_loci/05_cyp6_consensus.R).
#   pairwise_nt_identity_from_alignment.py <alignment_dir> > pairwise_gene_nt_identity.tsv
from pathlib import Path
import itertools
import sys

def read_fasta(fp):
    header = None
    seq = []
    with open(fp) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)

def pairwise_identity(seq1, seq2):
    matches = 0
    compared = 0
    for a, b in zip(seq1, seq2):
        if a == "-" or b == "-":
            continue
        if a == "N" or b == "N":
            continue
        compared += 1
        if a == b:
            matches += 1
    if compared == 0:
        return 0.0, 0
    return 100.0 * matches / compared, compared

def simplify_name(header):
    return header.split(".", 1)[0]

aln_dir = Path(sys.argv[1])

print("gene\tseq1\tseq2\tpct_identity\tcompared_sites")
for fp in sorted(aln_dir.glob("*.gene.aln.fa")):
    seqs = list(read_fasta(fp))
    gene = fp.name.replace(".gene.aln.fa", "")
    for (h1, s1), (h2, s2) in itertools.combinations(seqs, 2):
        pid, comp = pairwise_identity(s1, s2)
        print(f"{gene}\t{simplify_name(h1)}\t{simplify_name(h2)}\t{pid:.3f}\t{comp}")
