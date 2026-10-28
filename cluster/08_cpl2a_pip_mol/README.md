# PIP–MOL divergence at CpL2a beyond SNPs (cluster step 11)

PoolSNP calls SNPs only, so PIP–MOL differences in indels and in sequence that
fails to map are invisible in the allele-frequency data. These two scripts look
at the three reference pools (ME1/PIP, FL1/QUI, CAL1/MOL) across chr2:11–16 Mb
and the CpL2a interval chr2:13.5–14.0 Mb (Supplementary Methods, "Indel calling
and per-base depth profiling at CpL2a"; Supplementary Results, "PIP–MOL
sequence divergence at CpL2a beyond SNPs").

| Script | What it does | Outputs (under `cpl2a_analysis/`) |
|---|---|---|
| `01_cpl2a_depth_mapq.sh` (~15 min) | per-base depth with `samtools depth -a -q 10 -Q 20`, MAPQ distribution, flagstat, zero-depth runs > 100 bp | `depth/<pool>_cpl2a_{wide,candidate}_depth.tsv`, `_mapq_dist.tsv`, `_flagstat.txt`, `_zero_runs.tsv` |
| `02_cpl2a_indels_sv.sh` (~30 min) | A: `bcftools mpileup -q 10 -Q 20 -d 500 -a FORMAT/AD,FORMAT/DP` + `bcftools call -mv --ploidy 1` jointly on the three BAMs (haploid model for near-inbred colonies), indels in the interval, PIP–MOL fixed-different indels (one GT 0, the other 1); B: chr2 calls from the existing DELLY BCFs; C: counts of paired reads with mates on other chromosomes or inserts > 10 kb | `indels/cpl2a_indels_all.vcf.gz`, `cpl2a_indels_candidate.tsv`, `cpl2a_indels_pip_mol_fixed.tsv`; `delly/<pool>_inversions_{chr2_all,cpl2a}.tsv`; `discordant/<pool>_cpl2a_discordant.bam` |

Both activate the `cpl2a` conda environment (`envs/cpl2a.yml`: bcftools 1.9,
samtools 1.9) and then try `module load samtools` and `module load bcftools`;
a module that loads (samtools 1.16, bcftools 1.13 on this cluster) comes first
on `PATH`. The bcftools build that called the indels is recorded in the header
of `cpl2a_indels_all.vcf.gz` (`##bcftools_callVersion`). The two scripts can
run in parallel. Submit from the project root with `PROJECT_ROOT` set (see
`cluster/config.sh`):

```bash
sbatch scripts/08_cpl2a_pip_mol/01_cpl2a_depth_mapq.sh
sbatch scripts/08_cpl2a_pip_mol/02_cpl2a_indels_sv.sh
```

Inputs: the filtered BAMs `samples/<pool>/samfile/pip_1021/<pool>.filtered.bam`
(`01_alignment/filter_index_flagstat.sh`), the pipe-form reference
(`REF_FASTA`), and `discordant_pair_analysis/delly/<pool>_inversions.bcf`
(`07_inversion/run_delly_all3.sh`).
