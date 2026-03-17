# glmmDMR Tutorial (Detailed Reference)

このドキュメントは、glmmDMR の各スクリプトで実行している処理内容、主要オプション、入出力ファイル形式を実務向けに整理した詳細リファレンスです。

対象スクリプト:
- `summarize_extractor.py`
- `BinomTest.py`
- `prepare_matrix.sh`
- `run_glmmDMR.R`
- `DMR_merge.R`
- `make_binned_methylation_bigwig.R` (optional)
- `make_binned_variance_bigwig.R` (optional; Python script)

## 1. Pipeline Summary

標準フロー:
1. 上流で Bismark による methylation extraction を実行
2. `summarize_extractor.py` で site 集計
3. `BinomTest.py` で site フィルタ
4. `prepare_matrix.sh` で window matrix 生成
5. `run_glmmDMR.R` で window-level GLMM 推定
6. `DMR_merge.R` で DMR 統合
7. 必要に応じて bigWig 可視化

## 2. Upstream analysis example (before glmmDMR)

以下は、`glmmDMR` に入力する前段でよく使う Bismark ベースの最小例です。

```bash
bismark \
  --bowtie2 \
  -p ${core} \
  -o ${out} \
  ${ref} \
  -1 ${read}/${fa}_1.trimed.fq.gz \
  -2 ${read}/${fa}_2.trimed.fq.gz

samtools view -@ ${core} -q 42 -b ${out}/${fa}_1.trimed_bismark_bt2_pe.bam | \
  samtools sort -n -@ ${core} -o ${out}/${fa}.Q42.bam

mkdir -p ${call}/${fa}

# PCR deduplication
deduplicate_bismark \
  --paired \
  --output_dir ${out} \
  --bam \
  ${out}/${fa}.Q42.bam

# Call methylC
bismark_methylation_extractor \
  --paired-end \
  --output ${call}/${fa} \
  --parallel ${core} \
  --buffer_size 30% \
  --gzip \
  --genome_folder ${ref} \
  ${out}/${fa}.Q42.deduplicated.bam
```

この出力（`CpG_*.txt.gz`, `CHG_*.txt.gz`, `CHH_*.txt.gz`）を、次節の `summarize_extractor.py` に入力します。


## 3. Script-by-Script Details

## 3.1 `summarize_extractor.py`

Purpose:
- Bismark extractor 出力 (`*.txt.gz`) を context 横断で統合し、site 単位の集計表へ変換。

Core processing:
1. 入力ディレクトリ内の `CpG_*.txt.gz`, `CHG_*.txt.gz`, `CHH_*.txt.gz` を走査。
2. ファイル名 suffix (`OT/OB/CTOT/CTOB`) から strand を再定義。
3. Bismark ヘッダ行を除外して読み込み。
4. `code` から `meth/unmeth` を判定して site ごとに合算。
5. context ごとに一時 TSV を作成後、全 context を結合して最終 `tsv.gz` を出力。

Options:
- `-i, --input` (required): extractor 出力ディレクトリ
- `-o, --output` (required): 出力 `*.tsv.gz`
- `--threads` (default: 4): 並列処理スレッド数
- `--keep-temp`: context 別一時ファイルを保持

Input format:
- 各行に Bismark extractor 5列相当を含む `*.txt.gz`
- 想定列: read_id, strand, chr, pos, code

Output format (`*_summarized_output.tsv.gz`):
- `chr` (string)
- `pos` (int, 1-based)
- `strand` (`+`/`-`)
- `meth` (int)
- `unmeth` (int)
- `context` (`CpG`/`CHG`/`CHH`)

Example:
```bash
python summarize_extractor.py \
  -i calls/sample \
  -o matrix/sample_summarized_output.tsv.gz \
  --threads 4
```

## 3.2 `BinomTest.py`

Purpose:
- site 単位で binomial test + BH-FDR 補正を行い、低信頼シグナルを抑制。

Core processing:
1. 入力 TSV を読み込み。
2. `coverage = meth + unmeth` を計算。
3. `--min_coverage` で site フィルタ。
4. 帰無確率 `null_prob` を決定。
5. 各 site に binomial test を実行。
6. 全 p 値に対して BH 法で FDR 補正。
7. 非有意 (`FDR > threshold`) の `meth` を 0 に置換。
8. 最終出力はコア列のみ保存。

Options:
- `-i, --input` (required): summarize_extractor 出力
- `-o, --output` (required): 出力 `*.tsv.gz`
- `--null_prob` (default: None): 帰無確率
- `--fdr_threshold` (default: 0.05): FDR 閾値
- `--min_coverage` (default: 0): 最低 coverage
- `--nonconv_chr` (default: None): 非変換率推定に使う染色体
- `--threads` (default: 4): 並列数

Input format:
- 必須列: `chr, pos, strand, meth, unmeth, context`

Output format (`*_binomtest_result.tsv.gz`):
- `chr, pos, strand, meth, unmeth, context`
- `pval/FDR` は内部計算に使用し、最終出力には含めない実装

Behavior note:
- 非有意 site でも行は保持し、`meth=0` にする仕様。

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

## 3.3 `prepare_matrix.sh`

Purpose:
- 2群比較用の sliding-window matrix を context 別に生成。

Core processing:
1. 引数を解析し、入力ファイル存在を確認。
2. FASTA/FAI から `bedtools makewindows` で窓を作成。
3. 各入力 TSV.gz を BED 化。
4. `bedtools intersect` で window と cytosine を突合。
5. context ごとに matrix を gzip 出力。

Options:
- `--fasta` (required): FASTA または FAI
- `--group1` (required): group1 の TSV.gz 群
- `--group2` (required): group2 の TSV.gz 群
- `--group_labels` (default: `group1 group2`): 群名
- `--window` (default: 300): window 幅
- `--slide` (default: 200): slide 幅
- `--output` (default: `./matrix_out`): 出力ディレクトリ
- `--tmpdir`: 一時ディレクトリ

Input format:
- BinomTest 出力 (`chr,pos,strand,meth,unmeth,context`)

Output format (`<g1>_<g2>_<ctx>_matrix.tsv.gz`):
- 列順は `bedtools intersect -wa -wb` 由来の固定13列
- 実質内容: window座標 + site座標/群/sample/strand/meth/unmeth/coverage

Example:
```bash
bash prepare_matrix.sh \
  --fasta /path/to/TAIR10.fasta \
  --group1 binom/WT_1_binomtest_result.tsv.gz binom/WT_2_binomtest_result.tsv.gz \
  --group2 binom/MT_1_binomtest_result.tsv.gz binom/MT_2_binomtest_result.tsv.gz \
  --group_labels WT MT \
  --window 500 \
  --slide 300 \
  --output prep_out
```

## 3.4 `run_glmmDMR.R`

Purpose:
- window 単位で GLMM を当てて p 値、delta、要約統計を推定。

Core processing:
1. matrix を読み込み、固定列セットへ正規化。
2. context, coverage, site数, replicate 数で段階フィルタ。
3. 必要に応じて `prefilter_delta` で高速事前除外。
4. `mode=aggregate/site` に応じて window データ構築。
5. `family=binom/beta` で alt/null モデルを `glmmTMB` 推定。
6. LRT による p 値、`delta = mean_rate(group2)-mean_rate(group1)` を算出。
7. バッチ並列で全 window を処理し TSV.gz 出力。

Options (major):
- `-i, --infile` (required): 入力 matrix
- `-o, --out_prefix` (required): 出力 prefix
- `-c, --context` (default: NULL): context フィルタ
- `--group1`, `--group2` (required): 比較群ラベル
- `--min_reps_g1`, `--min_reps_g2` (default: 2): 群ごとの最小 replicate
- `--family` (default: `binom`): `binom` または `beta`
- `--mode` (default: `aggregate`): `aggregate` or `site`
- `--random_effect` (default: TRUE): `(1|sample)` を使用
- `--min_cov` (default: 0): site coverage フィルタ
- `--min_sites_win` (default: 0): window 内最小 site 数
- `--prefilter_delta` (default: 0): 事前デルタフィルタ
- `--workers` (default: 4): 並列 worker
- `--batches` (default: 50): 分割バッチ数
- `--max_globals_mb` (default: 1000): future global size
- `--seed` (default: 1): 乱数 seed

Input format:
- 必須列: `chr,start,end,cytosine_chr,cytosine_start,cytosine_end,context,group,sample,strand,meth,unmeth,coverage`

Output format (`*_fit_<family>_<mode>.tsv.gz`):
- `chr`
- `start`
- `end`
- `model`
- `p`
- `delta`
- `mean_rate1`
- `mean_rate2`
- `aic_diff`
- `bic_diff`

Example:
```bash
Rscript run_glmmDMR.R \
  -i prep_out/WT_MT_CpG_matrix.tsv.gz \
  -o glmm_out/WT_MT_CpG \
  --group1 WT --group2 MT \
  --family beta \
  --mode aggregate \
  --min_reps_g1 2 --min_reps_g2 2 \
  --min_sites_win 1 \
  --min_cov 5 \
  --random_effect \
  --workers 8 --batches 200 --max_globals_mb 2000
```

## 3.5 `DMR_merge.R`

Purpose:
- window-level 有意シグナルを統合して DMR を構築。

Core processing:
1. GLMM 出力を読み込み、`direction` を付与。
2. 指定 merge-mode で DMR 候補を構築。
3. 必要に応じて edge trim / post-filter / overlap merge。
4. DMR ごとの delta 要約を追加。
5. mode 別 TSV と BED を出力。

Supported merge modes:
- `Simes`
- `Stouffer`
- `single_seed`
- `multi_seed`
- `hybrid_seed`

Options (major):
- `--windows` (required): GLMM window 結果
- `--out-prefix` (default: `results/dmr`): 出力 prefix
- `--merge-mode` (default: `hybrid_seed`): 上記5種のみ有効
- `--p-seed` (default: 0.05)
- `--p-extend` (default: 0.01)
- `--max-gap-bp` (default: 200)
- `--min-windows` (default: 2)
- `--merge-overlaps` / `--merge-overlaps-gap`
- `--post-filter`, `--min-median-p`, `--min-consistent-frac`
- `--min-delta`, `--adaptive-delta`, `--adaptive-delta-method`, `--adaptive-delta-ratio`
- `--max-p-degradation`, `--max-final-p`, `--min-strong-windows`
- `--trim-weak-edges`
- `--min-dmr-length`
- `--max-median-p`
- `--seed-min-windows` (multi_seed向け)

Input format:
- 必須列: `chr, start, end, p, delta`

Output format (mode別):
- TSV: `*_dmrs_<mode>.tsv`
- BED: `*_dmrs_<mode>.bed`
- 主な列: `chr,start,end,n_windows,direction,combined_p`
- mode により `strong_frac`, `simes_p/stouffer_p`, `delta_mean/delta_max` が付加

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

## 3.6 `make_binned_methylation_bigwig.R` (optional)

Purpose:
- site methylation から bin ごとの平均 methylation を bigWig 化。

Core processing:
1. 入力 TSV.gz を読み込み、context で抽出。
2. `rate = meth/(meth+unmeth)` を算出。
3. bin に割り当てて平均値を集約。
4. GRanges 化し、genome size と整合を取って bw 出力。

Options:
- `-i, --input` (required): binomtest 結果
- `-b, --binsize` (default: 50): bin 幅
- `-o, --output` (required): 出力 bigWig
- `--genome` (required): chrom.sizes
- `--context` (default: `CpG`): `CpG/CHG/CHH`

Input format:
- 必須列: `chr,pos,meth,unmeth,context`

Output format:
- bigWig スコア = bin 平均 methylation rate

Example:
```bash
Rscript make_binned_methylation_bigwig.R \
  -i binom/sample_binomtest_result.tsv.gz \
  -b 50 \
  -o wig/sample_CpG.bw \
  --genome /path/to/TAIR10.chromInfo \
  --context CpG
```

## 3.7 `make_binned_variance_bigwig.R` (optional)

Purpose:
- 複数 bigWig の bin 平均値から replicate 間分散を計算して bigWig 化。

Note:
- ファイル名は `.R` ですが、実装は Python です。`python` で実行してください。

Core processing:
1. 入力 bigWig 群を同一 header 前提で読み込み。
2. 各 chromosome を bin 化して track ごとの平均値配列を取得。
3. 必要なら `log2p1` 正規化。
4. bin ごとに分散を計算。
5. `min-tracks` 未満の bin は NA として除外。
6. 出力 bigWig へ書き込み。

Options:
- `--inputs` (required): 入力 bigWig 群
- `--output` (required): 出力 bigWig
- `--bin-size` (default: 200)
- `--min-tracks` (default: 2)
- `--norm` (default: `none`): `none` or `log2p1`

Output format:
- bigWig スコア = bin ごとの track 間分散

Example:
```bash
python make_binned_variance_bigwig.R \
  --inputs wig/rep1_CpG.bw wig/rep2_CpG.bw wig/rep3_CpG.bw wig/rep4_CpG.bw \
  --output wig/group_CpG.variance.bw \
  --bin-size 200 \
  --min-tracks 2 \
  --norm none
```


## 4. Troubleshooting

1) No windows after filtering
- 原因: `--min_cov`, `--min_sites_win`, replicate 条件が厳しすぎる
- 対応: 閾値を段階的に緩める

2) GLMM が遅い / メモリ不足
- 対応: `--workers` を下げる、`--batches` を増やす、`TMPDIR` を高速かつ容量十分な場所へ

3) bigWig 生成で chromosome mismatch
- 原因: 入力 chr 名と chrom.sizes の不一致
- 対応: naming を統一して再実行
