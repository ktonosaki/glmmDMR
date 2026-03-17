<img width="688" height="192" alt="image" src="https://github.com/user-attachments/assets/69e3cb66-97c3-4ea0-9b10-f512209b0fef" />

glmmDMR is a small pipeline collection for DNA methylation analysis with window-level GLMM and DMR integration.

This repository contains scripts for:
- methylation extraction summarization
- site-level binomial filtering
- window matrix construction
- GLMM fitting
- DMR integration (multiple merge strategies)
- bigWig generation (mean and variance)

## Repository scripts

- summarize_extractor.py
  Summarize Bismark methylation extractor output (*.txt.gz) to per-site counts.
- BinomTest.py
  Per-site binomial test + FDR correction.
- prepare_matrix.sh
  Build sliding-window matrix from two groups of samples.
- run_glmmDMR.R
  Fit GLMM per window (binomial/beta, aggregate/site mode).
- DMR_merge.R
  Integrate significant windows into DMRs with multiple strategies.
- make_binned_methylation_bigwig.R
  Generate binned methylation bigWig.
- make_binned_variance_bigwig.R
  Generate binned variance bigWig from multiple bigWigs.
- simulate_site_data.R
  Simulate site-level methylation data for benchmarking.

## Install glmmDMR

Detailed step-by-step guide is available in tutorial_glmmDMR.md.

### Clone repository

```bash
git clone https://github.com/ktonosaki/glmmDMR.git
cd glmmDMR
```

### Make scripts executable (optional)

```bash
chmod +x summarize_extractor.py BinomTest.py prepare_matrix.sh
chmod +x run_glmmDMR.R DMR_merge.R
chmod +x make_binned_methylation_bigwig.R make_binned_variance_bigwig.R
```

### Quick check

```bash
python summarize_extractor.py --help
python BinomTest.py --help
Rscript run_glmmDMR.R --help
Rscript DMR_merge.R --help
```

### Run without installation

You can run scripts directly from the cloned repository (no package build required):

```bash
git clone https://github.com/ktonosaki/glmmDMR.git
cd glmmDMR
python summarize_extractor.py --help
Rscript run_glmmDMR.R --help
```



## Dependencies

Core tools used across scripts:
- Python: pandas, numpy, scipy, statsmodels, pyBigWig, tqdm
- R: optparse, data.table, glmmTMB, future, future.apply, GenomicRanges, rtracklayer
- CLI tools: bedtools, samtools, zcat

Recommended additional tools for upstream processing:
- fastp (read QC / trimming)
- SRA Toolkit (fasterq-dump)
- Bismark (+ bowtie2) for alignment and methylation extraction
- deepTools (e.g. bigwigAverage in replicate summaries)

### Install Python packages

```bash
python -m pip install pandas numpy scipy statsmodels pyBigWig tqdm
```

### Install R packages

```r
install.packages(c("optparse", "data.table", "glmmTMB", "future", "future.apply"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer"))
```

### Install CLI tools (example)

```bash
# macOS (Homebrew)
brew install bedtools samtools gzip

# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y bedtools samtools gzip
```

### Temporary directory setting (recommended)

Large datasets can generate large temporary files. Set TMPDIR to a high-capacity location:

```bash
export TMPDIR=/path/to/large_storage/tmp
mkdir -p "$TMPDIR"
```

To persist this setting, add it to your shell startup file (for example `~/.bashrc` or `~/.zshrc`).

You may also need Bismark/Fastp/SRA Toolkit depending on your upstream preprocessing.

## Typical workflow

The full workflow is based on the processing order used in 02.publis_data.sh:

1. Download and preprocess reads
   (for example: fasterq-dump, fastp, alignment and methylation extraction by Bismark)
2. Summarize methylation calls per sample
3. Run per-site binomial test
4. Convert to binned bigWig (optional)
5. Build group-wise sliding-window matrix
6. Run GLMM (run_glmmDMR.R)
7. Merge windows into DMRs (DMR_merge.R)
8. Compare with other methods / annotate genomic contexts (optional)

### Step 1: Download and preprocess reads

Goal: Generate high-quality alignment-ready FASTQ and methylation call files.

What is usually done:
- Download FASTQ files (for example with SRA Toolkit `fasterq-dump`).
- Perform read QC/trimming (for example with `fastp`).
- Align reads to the reference genome (typically Bismark + bowtie2).
- Remove PCR duplicates and run methylation extractor.

Expected output:
- Sample-wise Bismark extractor outputs (`*.txt.gz`) for each context/strand.

Example (paired-end, minimal flow):

```bash
# 1) Download FASTQ from SRA
fasterq-dump SRRXXXXXXX --split-files --threads 8 --outdir raw

# 2) Compress FASTQ
pigz -p 8 raw/SRRXXXXXXX_1.fastq
pigz -p 8 raw/SRRXXXXXXX_2.fastq

# 3) Trim reads
fastp \
  -i raw/SRRXXXXXXX_1.fastq.gz \
  -I raw/SRRXXXXXXX_2.fastq.gz \
  -o clean/sample_1.trimmed.fq.gz \
  -O clean/sample_2.trimmed.fq.gz \
  --thread 8

# 4) Align with Bismark (bowtie2)
bismark --bowtie2 -p 8 -o align /path/to/bismark_genome \
  -1 clean/sample_1.trimmed.fq.gz \
  -2 clean/sample_2.trimmed.fq.gz

# 5) Deduplicate and extract methylation calls
deduplicate_bismark --paired --output_dir align --bam align/sample_bismark_bt2_pe.bam
bismark_methylation_extractor --paired-end --gzip --parallel 8 \
  --genome_folder /path/to/bismark_genome \
  --output calls/sample \
  align/sample_bismark_bt2_pe.deduplicated.bam
```

### Step 2: Summarize methylation calls per sample

Goal: Convert extractor outputs into a compact per-site methylation table.

Script:
- `summarize_extractor.py`

Input:
- One sample extractor directory containing context files.

Output:
- `*_summarized_output.tsv.gz` with columns:
  `chr, pos, strand, meth, unmeth, context`.

Why this step matters:
- Standardizes output format for downstream filtering and matrix building.

### Step 3: Run per-site binomial test

Goal: Reduce likely non-informative methylation signal before window-level modeling.

Script:
- `BinomTest.py`

Input:
- `*_summarized_output.tsv.gz`.

Output:
- `*_binomtest_result.tsv.gz` with the same core columns.

Behavior note:
- Sites above the FDR threshold are retained but `meth` is set to 0.

### Step 4: Convert to binned bigWig (optional)

Goal: Create genome-browser tracks for quick visualization/QC.

Scripts:
- `make_binned_methylation_bigwig.R` (mean methylation per bin)
- `make_binned_variance_bigwig.R` (replicate variance per bin)

Input:
- Binomial-filtered TSVs and/or per-sample bigWig files.

Output:
- Context-specific bigWig tracks (`*.bw`).

### Step 5: Build group-wise sliding-window matrix

Goal: Convert site-level tables into the window matrix used by GLMM.

Script:
- `prepare_matrix.sh`

Input:
- Group1 and Group2 sets of `*_binomtest_result.tsv.gz`.
- Reference FASTA or FAI for window generation.

Output:
- Context-specific matrix files (`*_matrix.tsv.gz`) for each comparison.

Key controls:
- `--window`, `--slide` (resolution/smoothing tradeoff).
- Group labels used by downstream models.

### Step 6: Run GLMM (run_glmmDMR.R)

Goal: Fit window-level statistical models and estimate group differences.

Script:
- `run_glmmDMR.R`

Input:
- Matrix file from Step 5.

Output:
- Per-window fit table: `*_fit_<family>_<mode>.tsv.gz`.

Important options:
- `--family` (`binom` or `beta`)
- `--mode` (`aggregate` or `site`)
- replicate filters (`--min_reps_g1`, `--min_reps_g2`)
- parallelization (`--workers`, `--batches`)

### Step 7: Merge windows into DMRs (DMR_merge.R)

Goal: Integrate significant windows into biologically interpretable DMRs.

Script:
- `DMR_merge.R`

Input:
- GLMM fit output with required columns `chr,start,end,p,delta`.

Output:
- DMR tables and BED files (depending on merge mode/options).

Common strategy parameters:
- `--merge-mode` (for example `single_seed`, `stouffer_multi_seed`, `hybrid_seed`)
- `--p-seed`, `--p-extend`, `--min-windows`
- optional post-filtering / overlap merge flags.

### Step 8: Compare with other methods / genomic annotation (optional)

Goal: Benchmark glmmDMR calls and interpret genomic context.

Typical analyses:
- Compare against DSS / methylKit / DMRfinder outputs.
- Intersect DMR BED files with promoter, gene body, TE, and intergenic annotations.
- Summarize overlap counts and context-specific enrichment patterns.

## Input and output format

### summarize_extractor.py
- Input: Bismark extractor directory containing `*.txt.gz`
- Output: gzipped TSV with columns `chr, pos, strand, meth, unmeth, context`

### BinomTest.py
- Input: summarized TSV.gz from summarize_extractor.py
- Output: filtered TSV.gz with the same core columns
- Behavior: non-significant sites are retained with `meth=0` for weighted-level downstream usage

### prepare_matrix.sh
- Input: two groups of `*_binomtest_result.tsv.gz` and reference FASTA/FAI
- Output: context-specific sliding-window matrix files in TSV.gz

### run_glmmDMR.R
- Input: window matrix TSV.gz
- Output: per-window fit table `*_fit_<family>_<mode>.tsv.gz`

### DMR_merge.R
- Input: GLMM result table with required columns `chr,start,end,p,delta`
- Output: merged DMR tables (and BED files depending on mode/options)

## Quick test

Run help commands first to verify runtime dependencies and script entry points:

```bash
python summarize_extractor.py --help
python BinomTest.py --help
bash prepare_matrix.sh --help
Rscript run_glmmDMR.R --help
Rscript DMR_merge.R --help
```

## Minimal examples

### 1) Summarize extractor output

```bash
python summarize_extractor.py \
  -i sample_methylation_extractor_dir \
  -o sample_summarized_output.tsv.gz \
  --threads 4
```

### 2) Binomial test

```bash
python BinomTest.py \
  -i sample_summarized_output.tsv.gz \
  -o sample_binomtest_result.tsv.gz \
  --nonconv_chr chloroplast \
  --min_coverage 5 \
  --fdr_threshold 0.05 \
  --threads 4
```

### 3) Prepare matrix for two groups

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

### 4) Run GLMM

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

### 5) Merge DMRs

```bash
Rscript DMR_merge.R \
  --windows glmm_out/WT_MT_CpG_fit_beta_aggregate.tsv.gz \
  --out-prefix dmr_out/WT_MT_CpG \
  --merge-mode stouffer_multi_seed \
  --p-seed 0.05 \
  --p-extend 0.05 \
  --min-windows 1
```

### 6) bigWig outputs

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
python make_binned_variance_bigwig.R \
  --inputs rep1_CpG.bw rep2_CpG.bw rep3_CpG.bw rep4_CpG.bw \
  --output group_CpG.variance.bw \
  --bin-size 200 \
  --min-tracks 2 \
  --norm none
```

## Notes

- The pipeline is optimized for two-group comparisons.
- Contexts (CpG, CHG, CHH) are handled independently.
- Large datasets can require substantial memory/CPU; tune workers, batches, and filtering thresholds.

## Credits

- Images generated by Google Gemini

## Citation

If you use this repository in your analysis, please cite the repository URL together with the method papers used in your workflow (GLMM/DMR integration and upstream tools such as Bismark).
