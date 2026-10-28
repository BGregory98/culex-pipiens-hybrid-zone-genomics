# Data for: Urban enrichment of *molestus* ancestry and opposing temperature clines in a metropolitan *Culex pipiens* hybrid zone

Benjamin P. Gregory, Arielle Arsenault-Benoit, Katherine L. Bell and Megan L. Fritz

This deposit holds the large derived inputs that the R analysis stage of the accompanying code repository reads: allele counts and frequencies, the filtered variant calls, BayPass output, the JHB2020 annotation lifted onto the *Cx. pipiens* idCulPipi1 assembly, and the individual-level ddRADseq genotypes. Smaller inputs and every output table and figure are in the code repository itself.

- Code: https://github.com/BGregory98/culex-pipiens-hybrid-zone-genomics
- Raw Pool-seq reads (27 pools): NCBI SRA, BioProject PRJNA1551493 (private until 1 Jun 2027 or publication)
- Reference assembly idCulPipi1 (Hesson et al.): ENA project PRJEB67967

## Use with the code

Unpack the deposit anywhere and point the code at it, for example in a `.Renviron` file at the root of the code repository:

```
CULEX_DATA=/path/to/this/deposit
```

Every script resolves deposit paths with `deposit_file("<path>")` (R/helpers/config.R) and stops if a file is missing. The layout below must be kept as is.

## Contents

Sizes are in bytes; MD5 checksums are in `MANIFEST.tsv`.

| Path | Size | Contents | Made by |
|---|---|---|---|
| `data/vcfs/v4/poolsnp.biallelic.nomiss.vcf.gz` (+ `.tbi`) | 1,067,731,774 | PoolSNP v1.0.1 calls on idCulPipi1 for the 27 pools (24 field, 3 reference): biallelic SNPs with no missing pool, 11,629,890 SNPs. FORMAT carries per-pool RD/AD. | `cluster/_toplevel/PoolSNP_v4.sh`, then `bcftools view -m2 -M2 -v snps` (cluster/README.md, step 4) |
| `data/allele_counts_v4.tsv` | 2,037,000,149 | The same SNPs as a table: chromosome, position, REF, ALT, then `RD:AD` for each pool in VCF column order. No header. | `cluster/04_baypass/vcf_to_allele_counts.sh` |
| `data/allele_frequencies_v4.csv` | 1,905,247,460 | Per-pool alternate-allele frequency (AD / (RD + AD), two decimals) for the same SNPs, with a cumulative genome position. | `R/01_filtering/01_allele_frequency_table.R` |
| `results/baypass/v4/baypass_core_envsubset_summary_pi_xtx.out` | 1,011,800,477 | BayPass 3.1 core model, per-SNP XtXst (rows in the SNP order of the frequency table). | `cluster/04_baypass/baypass_core_envsubset.sh` |
| `results/baypass/v4/filtered_baypass_data.csv` | 149,912,034 | SNPs that are at least minor outliers for XtXst or for a Bayes factor, with all predictor BFs and outlier classes. | `cluster/04_baypass/filter_baypass_results.R` |
| `data/quinq_on_pip.gff3_polished` | 55,766,241 | The annotation of the *Cx. quinquefasciatus* JHB2020 assembly (GCA_015732765.1) lifted onto idCulPipi1 with Liftoff (Shumate and Salzberg 2021); `|` in feature names replaced with `_`. | Liftoff (Supplementary Methods) |
| `data/reference_genome/cxpip_genome.sanitized.fa` (+ `.fai`) | 541,394,902 | Chromosomes 1-3 of idCulPipi1 with `|` in sequence names replaced by `_`. | `cluster/02_mpileup/sanitize_headers.sh` |
| `data/reference_genome/chr1_only/chr1_quinq_vs_pip_chr1.asm10.paf` | 906,142 | minimap2 2.30 (`-x asm10`) alignment of chromosome 1 of JHB2020 against chromosome 1 of idCulPipi1. | `cluster/08_phylogeny/chr1_minimap.sh` |
| `data/reference_genome/chr1_only/chr1_quinq.gff3`, `chr1_pipiens.gff3` | 11,078,088; 11,355,479 | Chromosome-1 gene annotations of the two assemblies for the microsynteny figure. | Subset of the published annotations |
| `data/ddradseq/gprob2c0_all.txt`, `gprob2c1_all.txt` | 33,346,627 each | Individual-level ddRADseq data (field collections and library preparation by A. Arsenault-Benoit; sequencing and genotype calling by K. L. Bell with the pipeline of Bell et al. 2024): posterior mean genotypes (0-2) at 7,258 SNPs for 573 individuals from entropy version 1 (Gompert et al. 2014; Shastry et al. 2021), K = 2, one file per MCMC chain. | K. L. Bell |
| `data/ddradseq/indiv_ids.txt` | 8,603 | Individual, population and collection-event codes in the row order of the genotype files. | K. L. Bell |
| `data/ddradseq/ddrad_sample_sheet.csv` | 56,703 | Library sample sheet: sample ID, barcode, site, PCR species call, urbanization class, collection date. | A. Arsenault-Benoit |
| `data/border_shapefiles/cb_2024_us_{state,county}_500k.*` | 4.6 MB; 16.7 MB (.shp) | US Census Bureau 2024 cartographic boundary files (1:500,000), map backgrounds for Figures 1 and S7. Public domain; https://www.census.gov/geographies/mapping-files/time-series/geo/cartographic-boundary.html | US Census Bureau |

## Not included

- The eight CHELSA V2.1 climatologies (1981-2010; bio04, bio05, bio06, bio07, bio11, bio18, fcf, gdd10; Karger et al. 2017), about 2.4 GB of global rasters, are read only by `R/00_setup/01_site_predictors.R`. Download them from https://chelsa-climate.org into `data/chelsa_layers/` of this deposit to rerun that script. Its outputs (the cropped stack `data/local_env_layers_2.tif` and the site predictor tables) are in the code repository, so no other script needs them.
- The curated chemosensory annotation of Menna et al. (2026) on idCulPipi1, which belongs to that publication. Only `R/08_loci/11_chemosensory_gene_reference.R` reads it (place it at `data/menna_on_pip.gff3_polished` to rerun that script); its output, the chemosensory gene reference, is in the code repository.
- Raw reads and alignments (SRA, above).

## Licence

The data in this deposit are dedicated to the public domain under CC0 1.0 (https://creativecommons.org/publicdomain/zero/1.0/). Please still cite the article and this deposit.

Provenance of files derived from public resources, so users can trace them:

- `data/reference_genome/cxpip_genome.sanitized.fa`: the *Cx. pipiens* idCulPipi1 assembly (Hesson et al.; ENA project PRJEB67967) with sequence names cleaned for the pipeline.
- `data/quinq_on_pip.gff3_polished`: the *Cx. quinquefasciatus* JHB2020 annotation (assembly GCA_015732765.1) lifted onto idCulPipi1 with Liftoff (Shumate and Salzberg 2021). The `|` characters in some feature names were replaced with `_` because they broke the code.
- `data/reference_genome/chr1_only/`: chromosome 1 subsets of the two annotations and a minimap2 alignment between the two chromosomes.
- `data/border_shapefiles/`: US Census Bureau cartographic boundary files (public domain).
