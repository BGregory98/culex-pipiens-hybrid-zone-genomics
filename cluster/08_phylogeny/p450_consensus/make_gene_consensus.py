#!/usr/bin/env python3
# Per-pool whole-gene consensus sequences of the six CpL1a CYP6 paralogs: for
# every BAM in bamlist.txt (the 33 sequenced pools) and every gene in
# focal6.gene_coords.tsv, bcftools mpileup | call -c, calls with DP >= 5 and
# QUAL >= 20, positions with depth < 5 masked as N, bcftools consensus.
# Supplementary Methods "Coding consequences, phylogenetic placement,
# microsynteny and consensus sequences of the CpL1a CYP6 paralogs".
# Run from the project root. Output: p450_cluster6/08_consensus_cds/gene_nt/
# The coordinate table (focal6.gene_coords.tsv) is in the repository at data/consensus/p450s/.
from pathlib import Path
import subprocess
import sys
import re

REF = "reference_genomes/cxpip_genome/cxpip_genome.fasta"
BAMLIST = "bamlist.txt"
GENE_TABLE = "p450_cluster6/08_consensus_cds/focal6.gene_coords.tsv"
OUTDIR = Path("p450_cluster6/08_consensus_cds/gene_nt")
MINDP = 5
MINQUAL = 20

def run(cmd):
    print("+", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)

def read_fasta(path):
    header = None
    seq = []
    with open(path) as f:
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

OUTDIR.mkdir(parents=True, exist_ok=True)
tmp_ref = OUTDIR / "tmp_ref"
tmp_vcf = OUTDIR / "tmp_vcf"
tmp_mask = OUTDIR / "tmp_mask"
for d in [tmp_ref, tmp_vcf, tmp_mask]:
    d.mkdir(parents=True, exist_ok=True)

genes = {}
with open(GENE_TABLE) as f:
    for line in f:
        gene, chrom, start, end, strand = line.strip().split("\t")
        genes[gene] = (chrom, int(start), int(end), strand)

bams = []
with open(BAMLIST) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        bams.append(line)

if not Path(REF + ".fai").exists():
    run(["samtools", "faidx", REF])

for bam in bams:
    bam_path = Path(bam)
    if not bam_path.exists():
        print(f"WARNING: missing BAM {bam}", file=sys.stderr)
        continue

    sample = re.sub(r"\.bam$", "", bam_path.name)

    if not Path(str(bam_path) + ".bai").exists() and not Path(str(bam_path).replace(".bam", ".bai")).exists():
        run(["samtools", "index", str(bam_path)])

    for gene, (chrom, start, end, strand) in genes.items():
        region = f"{chrom}:{start}-{end}"
        tag = f"{sample}.{gene}"

        region_ref = tmp_ref / f"{tag}.ref.fa"
        vcf = tmp_vcf / f"{tag}.vcf.gz"
        mask = tmp_mask / f"{tag}.lowdp.bed"
        outfa = OUTDIR / f"{tag}.gene.fa"

        p = subprocess.run(
            ["samtools", "faidx", REF, region],
            check=True, capture_output=True, text=True
        )
        with open(region_ref, "w") as out:
            out.write(p.stdout)

        mp = subprocess.Popen(
            ["bcftools", "mpileup", "-f", REF, "-r", region, "-Ou", str(bam_path)],
            stdout=subprocess.PIPE
        )
        call = subprocess.Popen(
            ["bcftools", "call", "-c", "-Ou"],
            stdin=mp.stdout,
            stdout=subprocess.PIPE
        )
        subprocess.run(
            ["bcftools", "view", "-i", f"DP>={MINDP} && QUAL>={MINQUAL}", "-Oz", "-o", str(vcf)],
            stdin=call.stdout,
            check=True
        )
        mp.stdout.close()
        call.stdout.close()
        run(["tabix", "-f", "-p", "vcf", str(vcf)])

        depth = subprocess.run(
            ["samtools", "depth", "-r", region, "-a", str(bam_path)],
            check=True, capture_output=True, text=True
        )
        with open(mask, "w") as out:
            for line in depth.stdout.strip().splitlines():
                if not line:
                    continue
                c, pos, dp = line.split("\t")
                pos = int(pos)
                dp = int(dp)
                if dp < MINDP:
                    out.write(f"{c}\t{pos-1}\t{pos}\n")

        cons = subprocess.run(
            ["bcftools", "consensus", "-f", str(region_ref), "-m", str(mask), str(vcf)],
            check=True, capture_output=True, text=True
        )
        with open(outfa, "w") as out:
            seq = []
            for line in cons.stdout.splitlines():
                if line.startswith(">"):
                    out.write(f">{sample}|{gene}\n")
                else:
                    out.write(line + "\n")
