# glmmDMR Tutorial (Detailed Reference)

This document is a practical, detailed reference for glmmDMR, covering what each script does, key options, and input/output formats.

Target scripts:
- `summarize_extractor.py`
- `BinomTest.py`
- `prepare_matrix.sh`
- `run_glmmDMR.R`
- `merge_window.R`
- `make_binned_methylation_bigwig.R` (optional)
- `make_binned_variance_bigwig.py` (optional)

## 1. Pipeline Summary

Standard workflow:
1. Run methylation extraction with Bismark (upstream)
2. Summarize sites with `summarize_extractor.py`
3. Filter sites with `BinomTest.py`
4. Build a window matrix with `prepare_matrix.sh`
5. Estimate window-level GLMM statistics with `run_glmmDMR.R`
6. Integrate windows into DMRs with `DMR_merge.R`
7. Optionally generate bigWig tracks for visualization

Execution unit guide:
- Steps 2-3 (Sections 3.1-3.2) are run per sample.
- Step 4 (Section 3.3) integrates multiple sample outputs by group.

## 2. Upstream analysis example (before glmmDMR)

Below is a minimal Bismark-based upstream example commonly used before feeding data into glmmDMR.

```bash
bismark \
  --bowtie2 \
  -p 8 \
  -o align \
  /path/to/bismark_genome \
  -1 clean/sample_1.trimmed.fq.gz \
  -2 clean/sample_2.trimmed.fq.gz

samtools view -@ 8 -q 42 -b align/sample_1.trimmed_bismark_bt2_pe.bam | \
  samtools sort -n -@ 8 -o align/sample.Q42.bam

mkdir -p calls/sample

# PCR deduplication
deduplicate_bismark \
  --paired \
  --output_dir align \
  --bam \
  align/sample.Q42.bam

# Call methylC
bismark_methylation_extractor \
  --paired-end \
  --output calls/sample \
  --parallel 8 \
  --buffer_size 30% \
  --gzip \
  --genome_folder /path/to/bismark_genome \
  align/sample.Q42.deduplicated.bam
```

Use these outputs (`CpG_*.txt.gz`, `CHG_*.txt.gz`, `CHH_*.txt.gz`) as input to `summarize_extractor.py` in the next section.

## 3. Script-by-Script Details

## 3.1 `summarize_extractor.py`
- Scans `CpG_*.txt.gz`, `CHG_*.txt.gz`, and `CHH_*.txt.gz` in the input directory, aggregates counts per site, writes temporary TSV files per context, then merges all contexts into a final `tsv.gz` file.

Options:
- `-i, --input` (required): extractor output directory
- `-o, --output` (required): output `*.tsv.gz`
- `--threads` (default: 4): number of worker threads
- `--keep-temp`: keep context-specific temporary files

Input format:
- `*.txt.gz` files containing Bismark extractor-like 5-column records
- Expected columns: read_id, strand, chr, pos, code

Output format (`*_summarized_output.tsv.gz`):
| Column | Type / Values | Description |
| --- | --- | --- |
| `chr` | string | Chromosome name |
| `pos` | int (1-based) | Genomic position |
| `strand` | `+` / `-` | Strand |
| `meth` | int | Methylated count |
| `unmeth` | int | Unmethylated count |
| `context` | `CpG` / `CHG` / `CHH` | Cytosine context |

Example:
```bash
python summarize_extractor.py \
  -i calls/sample \
  -o matrix/sample_summarized_output.tsv.gz \
  --threads 4
```

## 3.2 `BinomTest.py`
Performs site-level binomial testing with BH-FDR correction to suppress low-confidence signals. Non-significant sites are zeroed.

Options:
- `-i, --input` (required): summarize_extractor output
- `-o, --output` (required): output `*.tsv.gz`
- `--nonconv_chr` (default: None): chromosome used to estimate non-conversion rate (recommended when available)
- `--null_prob` (default: None): null probability (use when `--nonconv_chr` cannot be used)
- `--fdr_threshold` (default: 0.05): FDR threshold
- `--min_coverage` (default: 0): minimum coverage
- `--threads` (default: 4): number of workers

Recommended usage for null probability:
- Preferred: set `--nonconv_chr` and estimate `null_prob` from a non-conversion control chromosome.
- Alternative: if no non-conversion control is available, provide a precomputed value via `--null_prob`.

Output format (`*_binomtest_result.tsv.gz`):
| Column | Type / Values | Description |
| --- | --- | --- |
| `chr` | string | Chromosome name |
| `pos` | int (1-based) | Genomic position |
| `strand` | `+` / `-` | Strand |
| `meth` | int | Methylated count (non-significant sites are set to 0) |
| `unmeth` | int | Unmethylated count |
| `context` | `CpG` / `CHG` / `CHH` | Cytosine context |

Example:
```bash
python BinomTest.py \
  -i matrix/sample_summarized_output.tsv.gz \
  -o binom/sample_binomtest_result.tsv.gz \
  --nonconv_chr chloroplast \
  --min_coverage 5 \
  --fdr_threshold 0.05 \
  --threads 4
```

Example (without non-conversion control chromosome; use precomputed null_prob):
```bash
python BinomTest.py \
  -i matrix/sample_summarized_output.tsv.gz \
  -o binom/sample_binomtest_result.tsv.gz \
  --null_prob 0.001 \
  --min_coverage 5 \
  --fdr_threshold 0.05 \
  --threads 4
```

## 3.3 `prepare_matrix.sh`

Purpose:
- Generates context-specific sliding-window matrices for two-group comparison.
- Takes per-sample `*_binomtest_result.tsv.gz` files from Sections 3.1-3.2 and integrates them via `--group1` and `--group2`.

Options:
- `--fasta` (required): FASTA or FAI
- `--group1` (required): group1 BinomTest TSV.gz files (multiple allowed)
- `--group2` (required): group2 BinomTest TSV.gz files (multiple allowed)
- `--group_labels` (default: `group1 group2`): group names
- `--window` (default: 300): window size
- `--slide` (default: 200): slide size
- `--output` (default: `./matrix_out`): output directory
- `--tmpdir`: temporary directory

Output format (`<g1>_<g2>_<ctx>_matrix.tsv.gz`):
| Item | Description |
| --- | --- |
| Layout | Fixed 13-column layout derived from `bedtools intersect -wa -wb` |
| Effective content | Window coordinates + site coordinates/group/sample/strand/meth/unmeth/coverage |

Example (2 samples per group):
```bash
bash prepare_matrix.sh \
  --fasta /path/to/TAIR10.fasta \
  --group1 binom/WT_1_binomtest_result.tsv.gz \
           binom/WT_2_binomtest_result.tsv.gz \
  --group2 binom/MT_1_binomtest_result.tsv.gz \
           binom/MT_2_binomtest_result.tsv.gz \
  --group_labels WT MT \
  --window 500 \
  --slide 300 \
  --output prep_out
```

Example (4 samples per group):
```bash
bash prepare_matrix.sh \
  --fasta /path/to/TAIR10.fasta \
  --group1 binom/WT_1_binomtest_result.tsv.gz \
           binom/WT_2_binomtest_result.tsv.gz \
           binom/WT_1_binomtest_result.tsv.gz \
           binom/WT_2_binomtest_result.tsv.gz \
  --group2 binom/MT_1_binomtest_result.tsv.gz \
           binom/MT_2_binomtest_result.tsv.gz \
           binom/MT_1_binomtest_result.tsv.gz \
           binom/MT_2_binomtest_result.tsv.gz \
  --group_labels WT MT \
  --window 500 \
  --slide 300 \
  --output prep_out
```

## 3.4 `run_glmmDMR.R`
- Fits GLMMs per window and estimates p-values, delta, and summary statistics. While multiple settings are available, `beta` family with `site` mode generally provides the highest accuracy. Stage-wise filtering by context, coverage, site count, and replicate count is supported. Optional `prefilter_delta` can accelerate computation by removing low-effect windows before fitting.

Options:

Input / Output:
- `-i, --infile` (required): input matrix
- `-o, --out_prefix` (required): output prefix

Data selection:
- `-c, --context` (default: NULL): context filter (`CpG`/`CHG`/`CHH`)
- `--group1`, `--group2` (required): comparison group labels

Model:
- `--family` (default: `beta`): `binom` or `beta` (`beta` + `site` is generally most accurate)
- `--mode` (default: `site`): `aggregate` (window-level summary) or `site` (site-level model)
- `--random_effect` (default: TRUE): use random effect `(1|sample)`

Pre-filter (before GLMM fitting):
- `--min_reps_g1`, `--min_reps_g2` (default: 2): minimum replicates per group
- `--min_cov` (default: 0): site coverage filter
- `--min_sites_win` (default: 0): minimum sites per window
- `--prefilter_delta` (default: 0): pre-remove windows with small delta (for speed)

Computation:
- `--workers` (default: 4): worker count
- `--batches` (default: 50): number of batches
- `--max_globals_mb` (default: 1000): future globals size limit (MB)
- `--seed` (default: 1): random seed

Output format (`*_fit_<family>_<mode>.tsv.gz`):
| Column | Description |
| --- | --- |
| `chr` | Chromosome name |
| `start` | Window start position |
| `end` | Window end position |
| `model` | Fitted model identifier |
| `p` | Statistical significance p-value |
| `delta` | Estimated methylation difference between groups |
| `mean_rate1` | Mean methylation rate in group1 |
| `mean_rate2` | Mean methylation rate in group2 |
| `aic_diff` | AIC difference metric |
| `bic_diff` | BIC difference metric |

Example:
```bash
Rscript run_glmmDMR.R \
  -i prep_out/WT_MT_CpG_matrix.tsv.gz \
  -o glmm_out/WT_MT_CpG \
  --group1 WT --group2 MT \
  --family beta \
  --mode site \
  --min_reps_g1 2 --min_reps_g2 2 \
  --min_sites_win 3 \
  --min_cov 5 \
  --random_effect \
  --workers 8 --batches 200 --max_globals_mb 2000
```

## 3.5 `merge_window.R`

- Integrates window-level significant signals into DMRs using the following three modes.

Supported merge modes:
- `single_seed`: starts from a strong single seed window and extends, conservative and easy to interpret.
- `multi_seed`: prioritizes linking regions containing multiple significant seed windows; strong for multi-peak regions.
- `hybrid_seed`: prioritizes `multi_seed` and complements uncovered regions with `single_seed`.

Options:

Input / Output:
- `--windows` (required): GLMM window results (`*_fit_<family>_<mode>.tsv.gz`)
- `--out-prefix` (default: `results/dmr`): output prefix

Mode:
- `--merge-mode` (default: `hybrid_seed`): only the three modes above are valid

Seed / Extension (core detection parameters):
- `--p-seed` (default: 0.05): p-value threshold for seed windows
- `--p-extend` (default: 0.05): p-value threshold for extension
- `--max-gap-bp` (default: 200): max gap to connect neighboring windows into one candidate
- `--min-windows` (default: 1): minimum windows required to report a DMR
- `--min-delta` (default: 0): minimum effect size threshold during extension
- `--max-p-degradation` (default: 1.2): allowed p-value worsening factor during extension (1.0 = no worsening)
- `--max-final-p` (default: 1.0): upper bound of final DMR combined p-value
- `--min-strong-windows` (default: 0.5): minimum fraction of windows with p <= p-seed in the final DMR

Adaptive delta threshold (automatic effect-size thresholding):
- `--adaptive-delta`: enable two-stage detection with data-driven delta thresholding
- `--adaptive-delta-method` (default: `median_ratio`): threshold method (`median_ratio`, `q50`, `q25`, `q10`, `mad`)
- `--adaptive-delta-ratio` (default: 0.6): ratio used in `median_ratio`

Notes:
- Using `adaptive-delta` in `multi_seed`, or in `hybrid_seed` (which runs `multi_seed` internally), often improves detection quality.
- In `method`, you can choose quantile-based thresholds (`q50`, `q25`, `q10`, etc.) based on data distribution.
- In `ratio`, you can directly control the threshold scaling.
- As an initial setting, `q25` (about 25%) is recommended.

Multi-seed specific:
- `--seed-min-windows` (default: 1): minimum seed windows for `multi_seed` and `hybrid_seed`

Post-filter (consistency checks after detection):
- `--post-filter`: enable quality filtering of candidate DMRs
- `--min-median-p` (default: 0.01): median p-value upper bound within a DMR (post-filter criterion)
- `--min-consistent-frac` (default: 0.5): minimum fraction of windows with p <= p-seed (post-filter criterion)

Length / median-p filters:
- `--min-dmr-length` (default: 0): minimum final DMR length (bp)
- `--max-median-p` (default: 1.0): independent filter on DMR median p-value

Overlap merge (post-detection re-merge):
- `--merge-overlaps`: re-merge overlapping/nearby same-direction DMRs
- `--merge-overlaps-gap` (default: 0): allowed gap between DMRs during re-merge

Output format (by mode):
| File | Pattern | Description |
| --- | --- | --- |
| TSV | `*_dmrs_<mode>.tsv` | DMR table output |
| BED | `*_dmrs_<mode>.bed` | Genome browser-compatible interval output |

Main TSV columns:

| Column | Description |
| --- | --- |
| `chr` | Chromosome name |
| `start` | DMR start position |
| `end` | DMR end position |
| `n_windows` | Number of windows included in the DMR |
| `direction` | Direction of methylation difference |
| `combined_p` | Combined p-value for the DMR |

Example:
```bash
Rscript merge_window.R \
  --windows glmm_out/WT_MT_CpG_fit_beta_aggregate.tsv.gz \
  --out-prefix dmr_out/WT_MT_CpG \
  --merge-mode hybrid_seed \
  --p-seed 0.05 \
  --p-extend 0.05 \
  --min-windows 1
```

## 3.6 `make_binned_methylation_bigwig.R` (optional)

- Converts site-level methylation (BinomTest TSV.gz output) into bin-wise mean methylation bigWig.

Options:
- `-i, --input` (required): BinomTest result
- `-b, --binsize` (default: 50): bin size
- `-o, --output` (required): output bigWig
- `--genome` (required): chrom.sizes
- `--context` (default: `CpG`): `CpG/CHG/CHH`

Example:
```bash
Rscript make_binned_methylation_bigwig.R \
  -i binom/sample_binomtest_result.tsv.gz \
  -b 50 \
  -o wig/sample_CpG.bw \
  --genome /path/to/TAIR10.chromInfo \
  --context CpG
```

## 3.7 `make_binned_variance_bigwig.py` (optional)
Computes replicate variance from bin-averaged values across multiple bigWig tracks and writes a variance bigWig. Optional normalization with `--norm log2p1` is available.

Options:
- `--inputs` (required): input bigWig files
- `--output` (required): output bigWig
- `--bin-size` (default: 200)
- `--min-tracks` (default: 2)
- `--norm` (default: `none`): `none` or `log2p1`

Output format:
| Item | Description |
| --- | --- |
| bigWig score | Per-bin variance across tracks |

Example:
```bash
python make_binned_variance_bigwig.py \
  --inputs wig/rep1_CpG.bw \
           wig/rep2_CpG.bw \
           wig/rep3_CpG.bw \
           wig/rep4_CpG.bw \
  --output wig/group_CpG.variance.bw \
  --bin-size 200 \
  --min-tracks 2 \
  --norm none
```

## 4. Troubleshooting

1) No windows after filtering
- Cause: thresholds such as `--min_cov`, `--min_sites_win`, or replicate conditions are too strict.
- Action: relax thresholds stepwise.

2) GLMM is slow / out of memory
- Action: reduce `--workers`, increase `--batches`, and set `TMPDIR` to a fast location with sufficient capacity.

3) Chromosome mismatch during bigWig generation
- Cause: chromosome names in inputs and chrom.sizes are inconsistent.
- Action: harmonize naming and rerun.
