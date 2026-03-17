# Tutorial: simulate_sites.R (CG-only)

This tutorial explains what `simulate_sites.R` does, how to run it, what each option means, and what files are produced.

The script simulates site-level methylation count data for two groups (for example, `WT` vs `MT`) with known DMR truth labels.

## 1. What This Script Generates

For a fixed context (`CG`), the script simulates:

1. Genomic DMR blocks with random directions (`hyper` or `hypo`) and effect sizes (`delta`).
2. Site positions across a synthetic chromosome.
3. Group- and replicate-specific methylation counts (`meth`, `unmeth`) with:
   - beta-binomial overdispersion,
   - replicate effects,
   - optional site-level heterogeneity,
   - missingness and coverage filtering.
4. Ground-truth labels for both blocks and sites.

## 2. Quick Start

Run with defaults:

```bash
Rscript simulate_sites.R
```

Run with custom settings:

```bash
Rscript simulate_sites.R \
  --out_dir results/site_window_sim \
  --chr_len 10000000 \
  --site_lambda 20000 \
  --block_frac 0.05 \
  --block_len_lo 300 \
  --block_len_hi 2000 \
  --delta_grid 0.1,0.2,0.3,0.4 \
  --groups WT,MT \
  --rep_per_group 4 \
  --base 0.7 \
  --dmr_site_prop 0.5 \
  --mu_site_sd 0.05 \
  --rho 0.3 \
  --logcov_mu 3.0 \
  --logcov_sd 1.2 \
  --miss_rate 0.30 \
  --min_cov 5 \
  --seed 202509
```

## 3. Model Summary

The script uses a CG-only simulation model.

- Base methylation level is set by `--base`.
- DMR effect (`delta`) is applied to group 2 only.
- If a site overlaps multiple DMR blocks, the delta with the largest absolute value is used.
- Per-site probability is sampled via beta distribution (with concentration derived from `rho`), then counts are drawn from binomial.

Informally:

```text
mu = base + delta * I(group == second_group) + replicate_effect + optional_site_effect
p ~ Beta(mu, concentration)
meth ~ Binomial(coverage, p)
unmeth = coverage - meth
```

## 4. Option Reference

## Required/General

`--out_dir`
- Output directory root.
- Default: `results/site_window_sim`

`--seed`
- Random seed for reproducibility.
- Default: `202509`

## Genome/Site Layout

`--chr_len`
- Simulated chromosome length in bp.
- Default: `5e6`

`--site_lambda`
- Site density per Mb for CG.
- Default: `4e4`

## DMR Block Construction

`--block_frac`
- Approximate fraction of chromosome covered by candidate DMR blocks.
- Default: `0.10`

`--block_len_lo`
- Minimum DMR block length (bp).
- Default: `500`

`--block_len_hi`
- Maximum DMR block length (bp).
- Default: `5000`

`--delta_grid`
- Comma-separated absolute effect sizes to sample from.
- Sign is assigned randomly (`hyper`/`hypo`).
- Default: `0.1,0.2,0.3`

## Group/Replicate Design

`--groups`
- Two group labels, comma-separated.
- Must contain exactly 2 groups.
- Default: `WT,MT`

`--rep_per_group`
- Number of replicates per group.
- Default: `3`

## Methylation Mean Structure

`--base`
- Baseline methylation rate for CG in group 1.
- Default: `0.7`

`--dmr_site_prop`
- Fraction of overlapping sites that receive DMR effect inside blocks.
- Range: `[0, 1]`
- Default: `1.0`

`--mu_site_sd`
- SD of optional site-level random effect added to mean methylation.
- Default: `0.0`

## Overdispersion and Coverage

`--rho`
- Base overdispersion scale for beta-binomial component.
- Default: `0.05`

`--logcov_mu`
- Mean of log-coverage distribution.
- Default: `log(20)`

`--logcov_sd`
- SD of log-coverage distribution.
- Default: `0.5`

`--miss_rate`
- Probability that a site is set to missing (`coverage = 0`) before coverage filtering.
- Default: `0.10`

`--min_cov`
- Minimum coverage threshold to keep a site.
- Default: `5`

## Context-specific Legacy Option

`--zeta_CHH`
- Zero-inflation rate for CHH.
- In the current script (`ctx = "CG"`), this option is not used in practice.
- Default: `0.40`

## 5. Outputs

All outputs are written to:

```text
<out_dir>/tsv/
```

with file names:

1. `sites_CG.tsv.gz`
2. `truth_blocks_CG.tsv.gz`
3. `truth_sites_CG.tsv.gz`

### 5.1 sites_CG.tsv.gz

Site-level simulated observations with truth labels joined.

Columns:
- `chr`: chromosome name (`chr1`)
- `pos`: genomic position (1-based)
- `sample`: sample ID (for example `WT01`)
- `group`: group label
- `replicate`: replicate index
- `context`: fixed as `CG`
- `meth`: methylated read count
- `unmeth`: unmethylated read count
- `truth`: `1` if site overlaps any truth block, else `0`
- `dir`: `hyper`, `hypo`, or `NA`

### 5.2 truth_blocks_CG.tsv.gz

Truth DMR blocks used during simulation.

Columns:
- `chr`
- `start`
- `end`
- `dir`
- `delta`

### 5.3 truth_sites_CG.tsv.gz

Unique sites with truth labels.

Columns:
- `chr`
- `pos`
- `truth`
- `dir`

## 6. Practical Notes

- The script is currently CG-only (`ctx <- "CG"`).
- `--groups` must provide exactly two labels.
- `--dmr_site_prop` must be in `[0, 1]`.
- For overlapping DMR blocks at one site, the script keeps the effect with largest absolute `delta`.

## 7. Suggested Parameter Sets

More realistic / moderate noise:

```bash
--site_lambda 20000 --block_frac 0.05 --rho 0.10 --miss_rate 0.10 --dmr_site_prop 0.8
```

Hard setting (more missingness/noise):

```bash
--site_lambda 20000 --block_frac 0.05 --rho 0.30 --miss_rate 0.30 --dmr_site_prop 0.5 --mu_site_sd 0.05
```

## 8. Minimal Validation After Run

```bash
zcat results/site_window_sim/tsv/sites_CG.tsv.gz | head
zcat results/site_window_sim/tsv/truth_blocks_CG.tsv.gz | head
zcat results/site_window_sim/tsv/truth_sites_CG.tsv.gz | head
```

