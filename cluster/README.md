# Cluster pipeline

The read-processing and genome-scan steps ran on a SLURM cluster. The scripts
here are published as they were run: every command, parameter and tool option
is the one that produced the files the R stage reads. They were not rerun for
publication. The R stage (`R/`) starts from their outputs, which are in the
data deposit or in this repository (see [Outputs read by the R stage](#outputs-read-by-the-r-stage)).

## Layout and configuration

On the cluster this directory was installed as `<project root>/scripts/`, and
jobs were submitted from the project root, so paths inside the scripts are
either relative to the project root (`locations.txt`, `logs/`,
`data/metadata/...`, `scripts/...`) or built from the variables in
`config.sh`. The scripts in `_toplevel/` lived in the project root itself.

`config.sh` holds the paths that were absolute in the scripts: the reference
FASTA files, the conda set-up and environment names, and the standalone
binaries (BWA-MEM2, Java/Picard, BayPass, DELLY, a samtools 1.21 build). It
derives all of them from one site-specific value, `PROJECT_ROOT`, which is not
stored in the repository. Set it once, either in the environment or in an
untracked `cluster/config.local.sh`:

```bash
# cluster/config.local.sh
PROJECT_ROOT=/path/to/project
```

With `PROJECT_ROOT` set to the original project directory, every path expands
to the one in the published run. Scripts source the file with
`source "${CLUSTER_CONFIG:-scripts/config.sh}"`; set `CLUSTER_CONFIG` to its
absolute path to submit from another directory.

Site-specific SLURM lines (`--mail-user`, `--mail-type`, `--account`,
`--partition`, `--qos`) were removed. Resource requests (`--time`, `--mem`,
`-c`/`--cpus-per-task`, `--array`) and log paths are kept as run; add your own
account, partition and QOS at submission if your cluster needs them.

Paths are as they were when each script ran. The per-pool directories
(`<pool>/samfile/pip_1021/...`) were later moved under `samples/`, which is
where the CpL2a scripts (step 11) find them.

Conda environments and pinned tool versions are in `../envs/` (see
[Tool versions](#tool-versions)).

## Run order

Pools: 24 field pools plus three reference pools, FL1 (QUI, *Cx.
quinquefasciatus*), ME1 (PIP, *Cx. pipiens* f. *pipiens*) and CAL1 (MOL, *Cx.
p.* f. *molestus*). Reference genome: *Cx. pipiens* idCulPipi1
(`cxpip_genome.fasta`, contig names in pipe form `ENA|OZ004311|OZ004311.1`;
`cxpip_genome.sanitized.fa`, the same sequence with `_` for `|`). Steps 1-2
ran on all 33 sequenced pools (`locations.txt`); from step 3 on, the 27 pools
in `locations_v4.txt` are used.

| Step | Scripts | What runs | Main outputs |
|---|---|---|---|
| 0 | (sequencing facility) | adapter and quality trimming, Trimmomatic v0.33 | `<pool>/ILLUMINA_DATA/*_trimmed.fastq.gz` |
| 1 Alignment | `01_alignment/combine_reads.sh`, `01_alignment/align_convert_sort.sh` | lane files concatenated; `bwa-mem2 mem -k 19 -t 8 -M -R <RG>`; samtools 1.21 sort, index, flagstat | `<pool>/samfile/pip_1021/<pool>.bam` |
| 2 Duplicates and filtering | `01_alignment/remove_duplicates.sh`, `01_alignment/filter_index_flagstat.sh` | Picard MarkDuplicates `REMOVE_DUPLICATES=true`; `samtools view -bh -f 2 -q 10` | `<pool>.dedup.bam` (DELLY input), `<pool>.filtered.bam` |
| 3 Pileup | `_toplevel/submit_v4.sh` → `_toplevel/mpileup_v4_array.sh`, `_toplevel/mpileup_v4_merge.sh` | per pool `samtools mpileup -B -q 10 -Q 20 -d 500 -a`, contig names sanitized inline; 27 pools pasted side by side | `pools27_v4.sanitized.mpileup.gz` |
| 4 Variant calling | `02_mpileup/run_sanitize.sh` → `02_mpileup/sanitize_headers.sh` (reference, once); `_toplevel/PoolSNP_v4.sh` (third job of `submit_v4.sh`) | PoolSNP `min-cov=10 max-cov=0.99 min-count=20 min-freq=0.05 miss-frac=0.2 base-quality=20 badsites=1 allsites=0`; then `bcftools view -e 'F_MISSING > 0'` | `PoolSNP_out_v4/poolsnp.vcf.gz`, `poolsnp.complete.vcf.gz` |
| 5 Biallelic filter | commands below | `bcftools view -m2 -M2 -v snps` | `PoolSNP_out_v4/poolsnp.biallelic.nomiss.vcf.gz` (11,629,890 SNPs) |
| 6 BayPass | `04_baypass/`: `vcf_to_allele_counts.sh`, `prepare_baypass_input_v4.R`, `baypass_core_envsubset.sh`, `baypass_env_IS.sh`, `pod_01_generate.R`, `pod_02_run_baypass.sh`, `pod_03_apply_threshold.R`, `run_baypass_filter.sh` → `filter_baypass_results.R` | core model and IS covariate model on the 24 field pools (`-nval 1000 -thin 25 -burnin 5000`, default 10 × 100 pilot runs); 100,000-SNP POD (seed 42) for the XtXst thresholds; join and coarse filter | `filtered_baypass_data.csv` |
| 7 Windowed f3 | `05_windowed_stats/`: `prepare_pooldata.slurm` → `prepare_pooldata.R`, `run_windowed_f3.sh` → `windowed_f3.R`, `combine_windowed_outputs.R` | poolfstat f3* in 1 Mb windows with ≥ 1,000 SNPs; BH within each reference pairing | `windowed_F3_all_comparisons_v4_1M.csv` |
| 8 Callable depth | `06_depth/glmm/` (see its README) | two passes over the pileup (99th-percentile cutoffs, then callable positions), CDS intersect, per-gene sums | `data/glmm_v4/*.tsv` |
| 9 Genome-wide null | `07_multilocus/genome_wide_null.slurm` → `01_genome_wide_null.R` | 10,000 width-matched random regions per locus, seed 42 | `genome_wide_null.rds` |
| 10 Inversion breakpoints | `_toplevel/insert_size_check.sh`, `07_inversion/run_delly_all3.sh` | insert-size summary; DELLY `call -t INV -q 20` on the reference pools' dedup BAMs, PASS calls with PE ≥ 5, INFO/MAPQ (median of the paired ends) ≥ 30, 50 kb–20 Mb | `discordant_pair_analysis/delly/<pool>_inversions.*.tsv` |
| 11 CpL2a PIP–MOL divergence | `08_cpl2a_pip_mol/` (see its README) | per-base depth and zero-depth runs; haploid indel calling on the three reference BAMs | `cpl2a_analysis/` |
| 12 CYP6 phylogeny, consensus, microsynteny | `08_phylogeny/run_p450_phylo.sh`, `08_phylogeny/p450_consensus/*.py`, `08_phylogeny/chr1_minimap.sh` | MAFFT L-INS-i/G-INS-i ± trimAl `-automated1`, IQ-TREE `-st AA -m MFP -alrt 1000 -B 1000`; per-pool CYP6 gene and CDS consensus with bcftools (`make_gene_consensus.py`, `make_cds_consensus.py`), then per-gene MAFFT alignment and pairwise identity; minimap2 `-x asm10` chr1 vs chr1 | `p450_cluster6/05_trees/all26.*`, `consensus/*.fa`, `chr1_quinq_vs_pip_chr1.asm10.paf` |

Steps 3–4 were submitted as one chain (`bash submit_v4.sh`). Steps 6, 7 and 9
need only step 5, and step 8 only step 3; they ran in parallel.

### Commands run interactively

All from the project root. Step 5 (the `bcftools view` command is recorded in
the VCF header):

```bash
module load bcftools            # bcftools 1.13
bcftools view -m2 -M2 -v snps -Oz \
  -o PoolSNP_out_v4/poolsnp.biallelic.nomiss.vcf.gz \
     PoolSNP_out_v4/poolsnp.complete.vcf.gz
bcftools index -t PoolSNP_out_v4/poolsnp.biallelic.nomiss.vcf.gz
```

Step 6 (the POD output prefix is the default, `baypass_pod`, as the POD run log
records):

```bash
bash scripts/04_baypass/vcf_to_allele_counts.sh \
  PoolSNP_out_v4/poolsnp.biallelic.nomiss.vcf.gz allele_counts_v4.tsv
module load R/4.5.3
Rscript scripts/04_baypass/prepare_baypass_input_v4.R \
  --counts allele_counts_v4.tsv --samples data/metadata/locations_v4.txt \
  --categories data/metadata/sample_categories.csv \
  --env data/metadata/site_coords_with_spatial_data.csv \
  --outdir SequencingData/data/baypass_v4
CORE=$(sbatch --parsable scripts/04_baypass/baypass_core_envsubset.sh)
sbatch --dependency=afterok:$CORE scripts/04_baypass/baypass_env_IS.sh
# after the core run:
Rscript scripts/04_baypass/pod_01_generate.R \
  --omega       SequencingData/data/baypass_v4/baypass_core_envsubset_mat_omega.out \
  --beta_params SequencingData/data/baypass_v4/baypass_core_envsubset_summary_beta_params.out \
  --counts      SequencingData/data/baypass_v4/baypass_STD_counts.txt \
  --poolsize    SequencingData/data/baypass_v4/baypass_poolsize_envsubset.txt \
  --baypass_utils software/baypass_public-master/utils/baypass_utils.R \
  --outdir      SequencingData/data/baypass_v4/pod
sbatch scripts/04_baypass/pod_02_run_baypass.sh SequencingData/data/baypass_v4/pod
Rscript scripts/04_baypass/pod_03_apply_threshold.R \
  --pod_xtx  SequencingData/data/baypass_v4/pod/baypass_pod_summary_pi_xtx.out \
  --real_xtx SequencingData/data/baypass_v4/baypass_core_envsubset_summary_pi_xtx.out \
  --outdir   SequencingData/data/baypass_v4
# after the IS run:
sbatch scripts/04_baypass/run_baypass_filter.sh SequencingData/data/baypass_v4
```

Step 7:

```bash
sbatch --export=WIN_SIZE=1000000 scripts/05_windowed_stats/prepare_pooldata.slurm
sbatch --array=1-24 --export=REF1=FL1,REF2=CAL1,WIN_SIZE=1000000 scripts/05_windowed_stats/run_windowed_f3.sh
sbatch --array=1-24 --export=REF1=ME1,REF2=CAL1,WIN_SIZE=1000000 scripts/05_windowed_stats/run_windowed_f3.sh
sbatch --array=1-24 --export=REF1=ME1,REF2=FL1,WIN_SIZE=1000000  scripts/05_windowed_stats/run_windowed_f3.sh
Rscript scripts/05_windowed_stats/combine_windowed_outputs.R 1000000
```

Steps 9–12 are submissions of the scripts named in the table
(`sbatch scripts/<dir>/<script>` from the project root), except the CYP6
alignments and trees of step 12, which were run on the cluster's submission
node rather than through SLURM, with 4 threads (host and thread count in the
MAFFT and IQ-TREE logs).

## Outputs read by the R stage

| Output on the cluster | Written by (step) | Read by | Location for the R stage |
|---|---|---|---|
| `<pool>.flagstat.txt`, `<pool>.filtered.flagstat.txt` → `mapping_rates_per_pool.csv` | 1, 2 (the table was compiled from the flagstat files; the command is not preserved) | `R/02_qc/01_mapping_rates.R` | repository `data/` |
| `PoolSNP_out_v4/poolsnp.biallelic.nomiss.vcf.gz` (+ `.tbi`) | 5 | `R/01_filtering/01_allele_frequency_table.R` and the VCF readers in `R/04_pop_structure/`, `R/07_multi_locus/`, `R/08_loci/` | deposit `data/vcfs/v4/` |
| `allele_counts_v4.tsv` | 6 (`vcf_to_allele_counts.sh`) | `R/07_multi_locus/06_allele_frequency_clines.R`, `10_windowed_pool_pca.R`, `24_background_ancestry_with_loci.R`; `R/08_loci/03_cyp6_coding_consequences.R`, `13_resistance_target_sites.R` | deposit `data/allele_counts_v4.tsv` |
| `SequencingData/data/baypass_v4/filtered_baypass_data.csv` | 6 (`filter_baypass_results.R`) | `R/03_baypass/01_baypass_outliers.R` | deposit `results/baypass/v4/` |
| `SequencingData/data/baypass_v4/baypass_core_envsubset_summary_pi_xtx.out` | 6 (core run) | `R/07_multi_locus/14_matched_null_xtx.R`, `15_matched_null_cpl2a_qui.R` | deposit `results/baypass/v4/` |
| `XtX_minor_T.rds`, `XtX_major_T.rds` (71.1249, 106.5979) | 6 (`pod_03_apply_threshold.R`) | values entered in `R/03_baypass/01_baypass_outliers.R` | — |
| `SequencingData/data/windowed_F3_all_comparisons_v4_1M.csv` | 7 | `R/07_multi_locus/08_windowed_f3.R` | repository `data/` |
| `data/glmm_v4/cds_depth_per_gene_pool.tsv`, `genome_callable_depth_per_pool.tsv` (and the `_all27` pair) | 8 | `R/06_depth/01_callable_depth_glmm.R`, `R/08_loci/07`–`09` | repository `data/glmm_v4/` |
| `genome_wide_null.rds` | 9 | `R/07_multi_locus/01_per_pool_ancestry.R`, `03_haplotype_block_diagnostics.R` | repository `results/multilocus/` |
| `discordant_pair_analysis/delly/CAL1_inversions.highconf.chr3_candidate.tsv` | 10 | CpL3a breakpoints, entered in `R/helpers/loci.R` | — |
| `p450_cluster6/05_trees/*.treefile` | 12 (`run_p450_phylo.sh`) | Figure S8 (`all26.linsi`, drawn in iTOL) | repository `data/consensus/p450s/trees/` |
| `p450_cluster6/08_consensus_cds/{gene_nt,cds_nt}/` → `pairwise_gene_nt_identity.tsv`, `pairwise_aa_identity.tsv` | 12 (`p450_consensus/make_*_consensus.py` with the gene and CDS coordinates in repository `data/consensus/p450s/focal6.*.tsv`, seqkit translation, per-gene MAFFT alignment, `p450_consensus/pairwise_*identity_from_alignment.py`) | `R/08_loci/05_cyp6_consensus.R` | repository `data/consensus/p450s/` |
| `reference_genomes/chr1_only/chr1_quinq_vs_pip_chr1.asm10.paf` | 12 (`chr1_minimap.sh`) | `R/08_loci/06_cyp6_microsynteny.R` | deposit `data/reference_genome/chr1_only/` |

Other numbers in the text come from the logs of these steps: the number of
biallelic SNPs (step 5), the callable-position count (155.1 M; step 8, pass 2
summary), the insert-size range (step 10) and the CpL2a indel and depth counts
(step 11).

The lifted annotation (`quinq_on_pip.gff3_polished`, Liftoff v1.6.3; deposit
`data/quinq_on_pip.gff3_polished`) is an input of step 8 and of the R stage.

## Tool versions

| Tool | Version | Source of the version | Used in step |
|---|---|---|---|
| Trimmomatic | 0.33 | sequencing facility | 0 |
| BWA-MEM2 | 2.2.1 | `bwa-mem2 version` | 1 |
| SAMtools (standalone build) | 1.21 | build directory name (`samtools-1.21`) | 1 |
| Picard (Java JDK 24.0.2) | 3.4.0 | `MarkDuplicates --version` | 2 |
| SAMtools / HTSlib (module) | 1.16.1 / 1.16 | default `module load samtools` ("Loading samtools/1.16" in the DELLY job log; VN:1.16.1 in the BAM @PG of the proper-pair filter) | 2, 3, 10, 12 |
| bcftools (module) | 1.13 | VCF headers (`1.13+htslib-1.16` in the v4 VCF; the CYP6 consensus scripts call `bcftools` from PATH, version not recorded; the CpL2a indel calls record 1.9 from the cpl2a env) | 4, 5, 12 |
| PoolSNP (conda `poolsnp`) | 1.0.1 (VCF header: PoolSnp-1.05) | `conda list -n poolsnp` | 4 |
| GNU parallel | 20170422 | `conda list -n poolsnp` | 4 |
| Python, numpy, tdigest (conda `poolsnp`) | 3.14.0, 1.26.4, 0.5.2.2 | `conda list -n poolsnp` | 8 |
| BayPass | 3.1 | `g_baypass -help` and run logs | 6 |
| R (module) | 4.5.3 | module name | 6, 7, 8, 9 |
| poolfstat, vcfR, data.table, dplyr, optparse | not recorded | — | 6, 7, 8, 9 |
| BEDTools (module) | not recorded | — | 8 |
| DELLY | 1.7.3 | `delly --version` | 10 |
| bcftools / SAMtools (conda `cpl2a`) | 1.9 / 1.9 | `conda list -n cpl2a` | 11 |
| MAFFT | 7.526 | alignment logs | 12 |
| trimAl | not recorded | — | 12 |
| IQ-TREE | 3.1.1 | tree logs | 12 |
| SeqKit | 2.13.0 | Supplementary Methods | 12 |
| minimap2 | 2.30 | `minimap2 --version` | 12 |

`../envs/poolsnp.yml` and `../envs/cpl2a.yml` reproduce the two conda
environments the scripts activate. The other files in `../envs/` pin the
versions above for the tools the scripts took from environment modules or
standalone builds; entries marked `# unconfirmed` are not recorded on the
cluster.
