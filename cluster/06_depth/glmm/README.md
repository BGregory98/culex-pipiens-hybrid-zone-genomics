# Callable-position read depth (cluster step 8)

Per-gene read depth at PoolSNP-callable positions, the input of the read-depth
GLMM (Methods, "Callable-position read-depth analysis"; Supplementary Methods,
same title). A position is callable when every one of the 27 pools has depth
≥ 10 (PoolSNP min-cov 10) and no pool exceeds its own per-chromosome 99th depth
percentile (PoolSNP max-cov 0.99).

## Scripts

| Step | Script (SLURM wrapper) | What it does | Output (`data/glmm_v4/`) |
|---|---|---|---|
| 1 | `01_compute_percentile_cutoffs.py` (`.slurm`) | pass 1 over the v4 pileup: per-pool, per-chromosome 99th-percentile depth (t-digest) | `pool_chrom_cutoffs_v4.tsv` |
| 2 | `02_extract_callable_depth.py` (`.slurm`) | pass 2: per-pool depth at every callable position | `callable_depth_v4.tsv.gz` |
| 3 | `03_intersect_cds.sh` (`.slurm`) | merged CDS per gene from the lifted annotation; callable positions inside CDS, tagged with gene_id | `cds_per_gene_merged.bed`, `callable_cds_depth.tsv.gz` |
| 4 | `04_aggregate_cds_depth.R` (`04_aggregate_cds_depth.slurm`) | per gene × field pool: summed depth and callable-position count (genes with ≥ 50 callable positions); genome-wide callable depth per pool | `cds_depth_per_gene_pool.tsv`, `genome_callable_depth_per_pool.tsv` |
| 4b | `04_aggregate_cds_depth.R --pool-set all` (`04b_aggregate_cds_depth_all27.slurm`) | step 4 with the three reference pools kept | `cds_depth_per_gene_pool_all27.tsv`, `genome_callable_depth_per_pool_all27.tsv` |

Inputs: `pools27_v4.sanitized.mpileup.gz` (`_toplevel/mpileup_v4_merge.sh`),
`data/metadata/locations_v4.txt`, `data/metadata/sample_categories.csv`,
`reference_genomes/quinq_on_pip.gff3_polished` (Liftoff annotation).

The four output tables are in the repository under `data/glmm_v4/` and are read
by `R/06_depth/01_callable_depth_glmm.R`, `R/08_loci/07_cyp6_copy_number.R`,
`R/08_loci/08_cyp6_depth_checks.R` and (`*_all27`)
`R/08_loci/09_cyp6_reference_pool_depth.R`.

## Running

From the project root, with `PROJECT_ROOT` set (see `cluster/config.sh`).
Steps 1 and 2 run in the `poolsnp` conda environment (`envs/poolsnp.yml`;
numpy and tdigest installed there with pip), step 3 uses the bedtools module
and step 4 the R/4.5.3 module with data.table and optparse.

```bash
mkdir -p data/glmm_v4 logs
J1=$(sbatch --parsable scripts/06_depth/glmm/01_compute_percentile_cutoffs.slurm)
J2=$(sbatch --parsable --dependency=afterok:$J1 scripts/06_depth/glmm/02_extract_callable_depth.slurm)
J3=$(sbatch --parsable --dependency=afterok:$J2 scripts/06_depth/glmm/03_intersect_cds.slurm)
J4=$(sbatch --parsable --dependency=afterok:$J3 scripts/06_depth/glmm/04_aggregate_cds_depth.slurm)
# reference-pool tables, after step 3 (check paths first with --dry-run):
bash  scripts/06_depth/glmm/04b_aggregate_cds_depth_all27.slurm --dry-run
sbatch scripts/06_depth/glmm/04b_aggregate_cds_depth_all27.slurm
```

Steps 1 and 2 each stream the whole pileup once (a few hours each); steps 3
and 4 take minutes. Step 3 caches `cds_per_gene_merged.bed` and reuses it if
present.

## Checks

- After step 2: every SNP in the v4 VCF should be a callable position.
- After step 3: the six CpL1a CYP6 paralogs and the CpL2a *Yan* ortholog
  (CQUJHB014704) should all appear in `cds_per_gene_merged.bed`.
- After step 4: the focal genes should pass the ≥ 50 callable-position filter.

## Notes

- Because the callable set requires depth ≥ 10 in all 27 pools, field and
  reference depths are measured at the same positions, and positions where a
  lineage drops out are excluded for every pool.
- `bedtools merge` in step 3 merges overlapping or abutting CDS of different
  genes into one interval labelled with both gene IDs (comma-separated), which
  step 4 treats as one unit.
