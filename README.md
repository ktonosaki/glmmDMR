(glmmDMR_logo.png)


## glmmDMR
glmmDMR partitions DNA methylation data into fixed-size sliding windows and applies a generalized linear mixed model (GLMM) with methylated cytosine counts per window as the response variable. By modeling the comparison target as a fixed effect while accounting for biological replicates, the framework estimates window-level methylation differences and statistical significance while incorporating between-replicate variability. To construct DMRs from significant windows, three integration modes are implemented: (i) `multi_seed`, which statistically combines multiple adjacent moderate signals, (ii) `single_seed`, which starts from a strongly significant single window and expands the region, and (iii) `hybrid_seed`, which prioritizes `multi_seed` and complements uncovered regions with `single_seed`.

Quick link:
- Detailed tutorial: [tutorial/tutorial_glmmDMR.md](tutorial/tutorial_glmmDMR.md)

## 1. Repository Contents

Core scripts:

- `summarize_extractor.py`: summarize Bismark extractor output (`*.txt.gz`) into per-site counts.
- `BinomTest.py`: per-site binomial test with FDR-based filtering behavior.
- `prepare_matrix.sh`: build two-group sliding-window matrices.
- `run_glmmDMR.R`: fit GLMM per window (`binom`/`beta`, `aggregate`/`site`).
- `DMR_merge.R`: merge significant windows into DMRs.
  - Supported `--merge-mode`: `single_seed`, `multi_seed`, `hybrid_seed`.
- `make_binned_methylation_bigwig.R` (optional): generate binned methylation bigWig.
- `make_binned_variance_bigwig.py` (optional): generate binned variance bigWig.

Simulation resources:

- `simulation/simulate_sites.R`
- `simulation/README.md`

Detailed method tutorial:

- `tutorial/tutorial_glmmDMR.md`

## 2. Installation

```bash
git clone https://github.com/ktonosaki/glmmDMR.git
cd glmmDMR
```

Optional executable flags:

```bash
chmod +x summarize_extractor.py BinomTest.py prepare_matrix.sh
chmod +x run_glmmDMR.R DMR_merge.R
chmod +x make_binned_methylation_bigwig.R make_binned_variance_bigwig.py
```

Python dependencies:

```bash
python -m pip install pandas numpy scipy statsmodels pyBigWig tqdm
```

R dependencies:

```r
install.packages(c("optparse", "data.table", "glmmTMB", "future", "future.apply"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer"))
```

Recommended CLI tools:

- bedtools
- samtools
- gzip / zcat
- Bismark (+ bowtie2) for upstream methylation extraction

## 3. Quick Environment Check

```bash
python summarize_extractor.py --help
python BinomTest.py --help
bash prepare_matrix.sh --help
Rscript run_glmmDMR.R --help
Rscript DMR_merge.R --help
```

## 4. Typical Workflow

1. Upstream methylation extraction (Bismark).
2. Summarize extractor output per sample.
3. Run per-site binomial filtering.
4. Build two-group sliding-window matrix.
5. Run GLMM per window.
6. Merge windows into DMRs.
7. Optionally generate bigWig tracks.

## 5. Step-by-Step Inputs and Outputs

### Step 1: `summarize_extractor.py`

Input:

- Bismark extractor output files (`*.txt.gz`) for one sample.

Output:

| File | Key columns |
| --- | --- |
| `*_summarized_output.tsv.gz` | `chr, pos, strand, meth, unmeth, context` |

Example:

```bash
python summarize_extractor.py \
  -i sample_methylation_extractor_dir \
  -o sample_summarized_output.tsv.gz \
  --threads 4
```

### Step 2: `BinomTest.py`

Input:

- `*_summarized_output.tsv.gz`

Output:

| File | Description |
| --- | --- |
| `*_binomtest_result.tsv.gz` | Site-level test result table (`chr, pos, strand, meth, unmeth, context`) |

Behavior note:

- Non-significant sites are retained, but `meth` is set to `0` for downstream usage.

Null probability note:

- Recommended: specify `--nonconv_chr` to estimate `null_prob` from a non-conversion control chromosome.
- If no non-conversion control chromosome is available: provide a precomputed value with `--null_prob`.

Example:

```bash
python BinomTest.py \
  -i sample_summarized_output.tsv.gz \
  -o sample_binomtest_result.tsv.gz \
  --nonconv_chr chloroplast \
  --min_coverage 5 \
  --fdr_threshold 0.05 \
  --threads 4
```

Alternative example (without non-conversion control chromosome):

```bash
python BinomTest.py \
  -i sample_summarized_output.tsv.gz \
  -o sample_binomtest_result.tsv.gz \
  --null_prob 0.012 \
  --min_coverage 5 \
  --fdr_threshold 0.05 \
  --threads 4
```

### Step 3: `prepare_matrix.sh`

Input:

- Group1 and Group2 sets of `*_binomtest_result.tsv.gz`
- Reference FASTA/FAI

Output:

| File pattern | Description |
| --- | --- |
| `<g1>_<g2>_<ctx>_matrix.tsv.gz` | Sliding-window matrix file |

Example:

```bash
bash prepare_matrix.sh \
  --fasta TAIR10.fasta \
  --group1 WT_1_binomtest_result.tsv.gz WT_2_binomtest_result.tsv.gz \
  --group2 MT_1_binomtest_result.tsv.gz MT_2_binomtest_result.tsv.gz \
  --group_labels WT MT \
  --window 500 \
  --slide 300 \
  --output prep_out
```

### Step 4: `run_glmmDMR.R`

Input:

- Window matrix TSV.GZ

Output:

| File pattern | Description |
| --- | --- |
| `*_fit_<family>_<mode>.tsv.gz` | Window-level GLMM result table |

Key options:

- `--family`: `binom` or `beta`
- `--mode`: `aggregate` or `site`
- `--workers`, `--batches`

Example:

```bash
Rscript run_glmmDMR.R \
  -i prep_out/WT_MT_CpG_matrix.tsv.gz \
  -o glmm_out/WT_MT_CpG \
  --group1 WT --group2 MT \
  --family beta \
  --mode site \
  --min_reps_g1 2 --min_reps_g2 2 \
  --min_sites_win 3 --min_cov 5 \
  --random_effect \
  --workers 8 --batches 200
```

### Step 5: `DMR_merge.R`

Input:

- GLMM fit table with required columns: `chr, start, end, p, delta`

Output:

| File pattern | Description |
| --- | --- |
| `*_dmrs_<mode>.tsv` | DMR table output |
| `*_dmrs_<mode>.bed` | BED interval output |

Supported merge modes:

- `single_seed`: start from a strong seed window and extend conservatively.
- `multi_seed`: prioritize regions that include multiple significant seed windows.
- `hybrid_seed`: run `multi_seed` first, then complement uncovered regions with `single_seed`.

Options by group:

Input / Output:

- `--windows` (required): input GLMM window result (`*_fit_<family>_<mode>.tsv.gz`)
- `--out-prefix` (default: `results/dmr`): output prefix

Mode:

- `--merge-mode` (default: `hybrid_seed`): `single_seed`, `multi_seed`, `hybrid_seed`

Seed / Extension:

- `--p-seed` (default: `0.05`): seed significance threshold
- `--p-extend` (default: `0.05`): extension threshold
- `--max-gap-bp` (default: `200`): maximum gap to connect adjacent windows
- `--min-windows` (default: `1`): minimum windows required for a DMR
- `--min-delta` (default: `0`): minimum effect size for extension
- `--max-p-degradation` (default: `1.2`): allowed p-value worsening during extension (`1.0` disables worsening)
- `--max-final-p` (default: `1.0`): maximum combined p-value of final DMR
- `--min-strong-windows` (default: `0.5`): minimum fraction of windows with `p <= p-seed`

Adaptive delta threshold:

- `--adaptive-delta`: enable automatic effect-size thresholding
- `--adaptive-delta-method` (default: `median_ratio`): `median_ratio`, `q50`, `q25`, `q10`, `mad`
- `--adaptive-delta-ratio` (default: `0.6`): ratio used for `median_ratio`

Adaptive delta usage note:

- `--adaptive-delta` is most useful in `multi_seed` or `hybrid_seed` workflows.
- Quantile methods (`q50`, `q25`, `q10`) can be selected based on data distribution.
- `q25` is a practical starting point for initial tuning.

Multi-seed specific:

- `--seed-min-windows` (default: `1`): minimum seed windows for `multi_seed` and `hybrid_seed`

Post-filter:

- `--post-filter`: enable post-detection DMR quality filtering
- `--min-median-p` (default: `0.01`): median p-value cutoff in post-filter
- `--min-consistent-frac` (default: `0.5`): minimum fraction with `p <= p-seed` in post-filter

Length / median-p filters:

- `--min-dmr-length` (default: `0`): minimum final DMR length (bp)
- `--max-median-p` (default: `1.0`): independent median p-value filter

Overlap merge:

- `--merge-overlaps`: re-merge overlapping/nearby DMRs with the same direction
- `--merge-overlaps-gap` (default: `0`): gap allowed for overlap re-merge

Example:

```bash
Rscript DMR_merge.R \
  --windows glmm_out/WT_MT_CpG_fit_beta_aggregate.tsv.gz \
  --out-prefix dmr_out/WT_MT_CpG \
  --merge-mode hybrid_seed \
  --p-seed 0.05 \
  --p-extend 0.05 \
  --min-windows 1
```

### Step 6: Optional bigWig tracks

Output:

| File | Description |
| --- | --- |
| `sample_CpG.bw` | Binned methylation bigWig |
| `group_CpG.variance.bw` | Binned variance bigWig |

Methylation bigWig:

```bash
Rscript make_binned_methylation_bigwig.R \
  -i sample_binomtest_result.tsv.gz \
  -b 50 \
  -o sample_CpG.bw \
  --genome TAIR10.chromInfo \
  --context CpG
```

Variance bigWig:

```bash
python make_binned_variance_bigwig.py \
  --inputs rep1_CpG.bw rep2_CpG.bw rep3_CpG.bw rep4_CpG.bw \
  --output group_CpG.variance.bw \
  --bin-size 200 \
  --min-tracks 2 \
  --norm none
```

## 6. Runtime Notes

- Pipeline assumes two-group comparisons.
- Contexts are usually processed independently.
- For large datasets, tune filtering and parallel options for memory/runtime balance.

For large jobs, set a high-capacity temporary directory:

```bash
export TMPDIR=/path/to/large_storage/tmp
mkdir -p "$TMPDIR"
```

## 7. Simulation and Benchmarking

For controlled benchmark data:

- simulation/README.md
- simulation/simulate_sites.R

## 8. Citation

If you use this repository, cite the repository URL and upstream method/tool papers used in your workflow.
