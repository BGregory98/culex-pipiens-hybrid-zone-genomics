#!/usr/bin/env python3
# Per-pool CDS consensus sequences of the six CpL1a CYP6 paralogs, built as in
# make_gene_consensus.py for each CDS piece in focal6.cds_coords.fixed.tsv,
# concatenated and reverse-complemented for minus-strand genes. The CDS
# consensus is then translated (seqkit) and aligned per gene (MAFFT
# --localpair --maxiterate 1000) for pairwise_identity_from_alignment.py.
# Run from the project root. Output: p450_cluster6/08_consensus_cds/cds_nt/
# The coordinate table (focal6.cds_coords.fixed.tsv) is in the repository at data/consensus/p450s/.
from pathlib import Path
import subprocess
import sys
import re

REF = "reference_genomes/cxpip_genome/cxpip_genome.fasta"
BAMLIST = "bamlist.txt"
CDS_TABLE = "p450_cluster6/08_consensus_cds/focal6.cds_coords.fixed.tsv"
OUTDIR = Path("p450_cluster6/08_consensus_cds")
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

def revcomp(seq):
    tbl = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tbl)[::-1]

for sub in ["tmp_ref", "tmp_vcf", "tmp_mask", "tmp_consensus", "cds_nt"]:
    (OUTDIR / sub).mkdir(parents=True, exist_ok=True)

cds = {}
strand = {}
with open(CDS_TABLE) as f:
    for line in f:
        gene, chrom, start, end, s = line.strip().split("\t")
        start = int(start)
        end = int(end)
        cds.setdefault(gene, []).append((chrom, start, end))
        strand[gene] = s

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

    for gene, pieces in cds.items():
        piece_seqs = []

        for i, (chrom, start, end) in enumerate(sorted(pieces, key=lambda x: x[1]), start=1):
            region = f"{chrom}:{start}-{end}"
            tag = f"{sample}.{gene}.cds{i}"

            region_ref = OUTDIR / "tmp_ref" / f"{tag}.ref.fa"
            vcf = OUTDIR / "tmp_vcf" / f"{tag}.vcf.gz"
            mask = OUTDIR / "tmp_mask" / f"{tag}.lowdp.bed"
            outfa = OUTDIR / "tmp_consensus" / f"{tag}.fa"

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
                out.write(cons.stdout)

            seq = None
            for h, s in read_fasta(outfa):
                seq = s
                break
            if seq is None:
                raise RuntimeError(f"No consensus sequence for {tag}")
            piece_seqs.append(seq)

        cds_seq = "".join(piece_seqs)
        if strand[gene] == "-":
            cds_seq = revcomp(cds_seq)

        final_out = OUTDIR / "cds_nt" / f"{sample}.{gene}.cds.fa"
        with open(final_out, "w") as out:
            out.write(f">{sample}|{gene}\n")
            for i in range(0, len(cds_seq), 80):
                out.write(cds_seq[i:i+80] + "\n")
