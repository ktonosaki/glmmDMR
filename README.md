
## glmmDMR
A sliding-window strategy was used to partition DNA methylation data into fixed-size windows, and a generalized linear mixed model (GLMM) was applied with methylated cytosine counts in each window as the response variable. By modeling the comparison target as a fixed effect while accounting for biological replicates, the framework estimates window-level methylation differences and statistical significance while incorporating between-replicate variability. In addition, to construct DMRs from windows judged significant, two integration methods were implemented: (i) a Multi-seed approach that statistically combines multiple adjacent moderate signals, and (ii) a Single-seed approach that starts from a strongly significant single window and expands the region.

## 1. Repository Contents

Core scripts:

- summarize_extractor.py: summarize Bismark extractor output (`*.txt.gz`) into per-site counts.
- BinomTest.py: per-site binomial test with FDR-based filtering behavior.
- prepare_matrix.sh: build two-group sliding-window matrices.
- run_glmmDMR.R: fit GLMM per window (`binom`/`beta`, `aggregate`/`site`).
- DMR_merge.R: merge significant windows into DMRs.
  - Supported `--merge-mode`: `Simes`, `Stouffer`, `single_seed`, `multi_seed`, `hybrid_seed`.
- make_binned_methylation_bigwig.R (optional): generate binned methylation bigWig.
- make_binned_variance_bigwig.R (optional): generate binned variance bigWig.

Simulation resources:

- simulation/simulate_sites.R
- simulation/README.md

Detailed method tutorial:

- tutorial/tutorial_glmmDMR.md

## 2. Installation

```bash
git clone https://github.com/ktonosaki/glmmDMR.git
cd glmmDMR
```

Optional executable flags:

```bash
chmod +x summarize_extractor.py BinomTest.py prepare_matrix.sh
chmod +x run_glmmDMR.R DMR_merge.R
chmod +x make_binned_methylation_bigwig.R make_binned_variance_bigwig.R
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

### Step 1: summarize_extractor.py

Input:

- Bismark extractor output files (`*.txt.gz`) for one sample.

Output:

- `*_summarized_output.tsv.gz`
- Key columns: `chr, pos, strand, meth, unmeth, context`

Example:

```bash
python summarize_extractor.py \
  -i sample_methylation_extractor_dir \
  -o sample_summarized_output.tsv.gz \
  --threads 4
```

### Step 2: BinomTest.py

Input:

- `*_summarized_output.tsv.gz`

Output:

- `*_binomtest_result.tsv.gz`

Behavior note:

- Non-significant sites are retained, but `meth` is set to `0` for downstream usage.

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

### Step 3: prepare_matrix.sh

Input:

- Group1 and Group2 sets of `*_binomtest_result.tsv.gz`
- Reference FASTA/FAI

Output:

- Sliding-window matrix TSV.GZ files

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

### Step 4: run_glmmDMR.R

Input:

- Window matrix TSV.GZ

Output:

- `*_fit_<family>_<mode>.tsv.gz`

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
  --mode aggregate \
  --min_reps_g1 2 --min_reps_g2 2 \
  --min_sites_win 1 --min_cov 5 \
  --random_effect \
  --workers 8 --batches 200
```

### Step 5: DMR_merge.R

Input:

- GLMM fit table with required columns: `chr, start, end, p, delta`

Output:

- DMR TSV files and BED files

Key options:

- `--merge-mode`: `Simes`, `Stouffer`, `single_seed`, `multi_seed`, `hybrid_seed`
- `--p-seed`, `--p-extend`, `--min-windows`

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
Rscript make_binned_variance_bigwig.R \
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
