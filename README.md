# Urban enrichment of *molestus* ancestry and opposing temperature clines in a metropolitan *Culex pipiens* hybrid zone

Code, small inputs, result tables and figures for

> Gregory, B. P., Arsenault-Benoit, A., Bell, K. L., & Fritz, M. L. Urban enrichment of *molestus* ancestry and opposing temperature clines in a metropolitan *Culex pipiens* hybrid zone. *Molecular Biology and Evolution* (submitted).

Whole-genome Pool-seq of 24 field pools from the Baltimore–Washington corridor and three reference pools (*Cx. quinquefasciatus*, QUI; *Cx. pipiens* f. *pipiens*, PIP; *Cx. pipiens* f. *molestus*, MOL) shows *molestus*-like ancestry rising with urbanization across the genome, and three environmentally structured loci with different donors, architectures and gradients:

| Locus | Interval (idCulPipi1) | Donor | Gradient |
|---|---|---|---|
| CpL1a | chr1 108.5–111.0 Mb (CYP6 cluster) | QUI ancestry rises with temperature | maximum summer temperature (bio05) |
| CpL2a | chr2 12.5–15.0 Mb block, 13.5–14.0 Mb core (*Yan*) | QUI ancestry excluded | bio05 |
| CpL3a | chr3 149,402,438–161,021,387 (inversion) | MOL ancestry rises with urbanization, in step with genome-wide MOL ancestry | percent impervious surface |

The coordinates used by every script are in `R/helpers/loci.R`.

## Layout

```
R/                 the analysis, rerunnable from the data deposit
  00_setup/        site predictors, Table S1, Figure 1, Figure S1
  01_filtering/    allele-frequency table from the VCF (optional)
  02_qc/           mapping-rate check
  03_baypass/      BayPass candidate SNPs, cluster inventory, Figure 3
  04_pop_structure/  FST, PCA, isolation by distance, ddRADseq PCA, diversity
  06_depth/        callable-depth GLMM, Figure 6
  07_multi_locus/  donor ancestry tests, haplotype blocks, LD, f3, windowed PCA, nulls
  08_loci/         locus-level analyses (CYP6 cluster, chemosensory genes, kdr)
  helpers/         configuration and shared definitions
  run_all.R        runs everything in order
cluster/           the HPC pipeline as it was run (reads -> VCF, BayPass, f3, depth, null, DELLY, phylogeny)
envs/              tool versions for the cluster pipeline
data/              small inputs, and tables written by the scripts
results/           result tables written by the scripts
figures/final/     every manuscript figure, written by the scripts (Figure S8 drawn in iTOL)
docs/              MANUSCRIPT_MAP.md (figure, table and number -> script), DATA_DEPOSIT.md
renv.lock          R package versions
```

## Rerunning the analysis

1. **R and packages.** R 4.6.1. From the repository root, `install.packages("renv"); renv::restore()` installs the package versions in `renv.lock` (CRAN and Bioconductor 3.23).
2. **Data.** Download the data deposit (https://doi.org/10.5061/dryad.gmsbcc367) and unpack it. It holds the allele counts and frequencies, the filtered VCF, BayPass output, lifted annotations, the reference sequence and the ddRADseq genotypes (about 6.5 GB; contents in `docs/DATA_DEPOSIT.md`). Tell the code where it is with a `.Renviron` file in the repository root (see `.Renviron.example`):
   ```
   CULEX_DATA=/path/to/culex-deposit
   ```
3. **Run.** From the repository root, start a new R session (so `.Renviron` is read) and run
   ```
   Rscript R/run_all.R
   ```
   or any single script, e.g. `Rscript R/07_multi_locus/01_per_pool_ancestry.R`. Scripts must be run from the repository root; each stops with a message if an input is missing. Run order is stage by stage (00, 01, 02, 03, 04, 06, 07, 08) and by file number within a stage; every script reads only outputs of earlier ones. The whole stage takes about 3 h on a 16 GB laptop; scripts that read the full allele-frequency table need about 8 GB of memory.

Every table in `results/` and `data/`, and every figure in `figures/final/`, is rewritten by the run. The repository versions are the ones reported in the paper, so `git diff` after a run shows any difference.

Three steps need more than the deposit, and their outputs are in the repository: `R/00_setup/01_site_predictors.R` extracts the site predictors from the global CHELSA climatologies (download instructions in `docs/DATA_DEPOSIT.md`); `R/08_loci/11_chemosensory_gene_reference.R` builds the chemosensory gene reference from the curated annotation of Menna et al. (2026), which belongs to that publication; `run_all.R` skips both when their inputs are absent. `R/01_filtering/01_allele_frequency_table.R` rebuilds the deposit's allele-frequency table from the VCF and is skipped by default.

## The cluster pipeline

`cluster/` holds the scripts that produced the deposit's inputs from raw reads, as run on a SLURM cluster: alignment (BWA-MEM2), pileup, PoolSNP variant calling, BayPass, windowed f3 (poolfstat), callable depth, the genome-wide null for the donor-ancestry test, DELLY inversion calls, and the CYP6 consensus sequences and phylogeny. `cluster/README.md` gives the run order, parameters and the outputs that feed the R stage; `envs/` gives tool versions. Paths are set in `cluster/config.sh`; site-specific SLURM account and partition lines were removed. These scripts are documentation of what was run and are not packaged to rerun elsewhere without adaptation.

## Data sources

Raw Pool-seq reads: NCBI SRA, BioProject PRJNA1551493. Reference assembly idCulPipi1: ENA PRJEB67967. Individual ddRADseq genotypes: generated by A. Arsenault-Benoit and K. L. Bell with the pipeline of Bell et al. (2024), in the data deposit. Climate: CHELSA V2.1 (Karger et al. 2017). Land cover and impervious surface: NLCD 2021. Map boundaries: US Census Bureau cartographic boundary files.

## Licence

Code: MIT (`LICENSE`). Data, result tables, figures and documentation: CC BY 4.0 (`LICENSE-data`). The data deposit on Dryad is released under CC0.

## Citation

See `CITATION.cff`.
