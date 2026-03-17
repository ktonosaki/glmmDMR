# glmmDMR Tutorial (Detailed Reference)

このドキュメントは、glmmDMR の各スクリプトで実行している処理内容、主要オプション、入出力ファイル形式を実務向けに整理した詳細リファレンスです。

対象スクリプト:
- `summarize_extractor.py`
- `BinomTest.py`
- `prepare_matrix.sh`
- `run_glmmDMR.R`
- `DMR_merge.R`
- `make_binned_methylation_bigwig.R` (optional)
- `make_binned_variance_bigwig.py` (optional)

## 1. Pipeline Summary

標準フロー:
1. 上流で Bismark による methylation extraction を実行
2. `summarize_extractor.py` で site 集計
3. `BinomTest.py` で site フィルタ
4. `prepare_matrix.sh` で window matrix 生成
5. `run_glmmDMR.R` で window-level GLMM 推定
6. `DMR_merge.R` で DMR 統合
7. 必要に応じて bigWig 可視化

実行単位の目安:
- ステップ 2〜3（= Section 3.1〜3.2）は「サンプルごと」に実行します。
- ステップ 4（= Section 3.3）で、複数サンプルの結果を群ごとにまとめて統合します。

## 2. Upstream analysis example (before glmmDMR)

以下は、`glmmDMR` に入力する前段でよく使う Bismark ベースの最小例です。

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

この出力（`CpG_*.txt.gz`, `CHG_*.txt.gz`, `CHH_*.txt.gz`）を、次節の `summarize_extractor.py` に入力します。


## 3. Script-by-Script Details

## 3.1 `summarize_extractor.py`
- 入力ディレクトリ内の `CpG_*.txt.gz`, `CHG_*.txt.gz`, `CHH_*.txt.gz` を走査し、site ごとに合算して、context ごとに一時 TSV を作成後、全 context を結合して最終 `tsv.gz` を出力。

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
site 単位で binomial test + BH-FDR 補正を行い、低信頼シグナルを抑制。非有意 site は0に置換する。

Options:
- `-i, --input` (required): summarize_extractor 出力
- `-o, --output` (required): 出力 `*.tsv.gz`
- `--nonconv_chr` (default: None): 非変換率推定に使う染色体（基本はこちらを推奨）
- `--null_prob` (default: None): 帰無確率（`--nonconv_chr` が使えない場合に、事前計算した値を指定）
- `--fdr_threshold` (default: 0.05): FDR 閾値
- `--min_coverage` (default: 0): 最低 coverage
- `--threads` (default: 4): 並列数

Recommended usage for null probability:
- 基本運用: `--nonconv_chr` を指定して、非変換コントロール染色体から `null_prob` を推定する。
- 代替運用: 非変換コントロールがない場合は、外部で計算した `null_prob` を `--null_prob` に明示指定して実行する。

Output format (`*_binomtest_result.tsv.gz`):
- `chr, pos, strand, meth, unmeth, context`

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

Example (without nonconv control chromosome; use precomputed null_prob):
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
- 2群比較用の sliding-window matrix を context 別に生成。
- Section 3.1〜3.2 で各サンプルごとに作成した `*_binomtest_result.tsv.gz` を、ここで `--group1` / `--group2` にまとめて渡して統合します。

Options:
- `--fasta` (required): FASTA または FAI
- `--group1` (required): group1 の BinomTest 出力 TSV.gz 群（複数指定）
- `--group2` (required): group2 の BinomTest 出力 TSV.gz 群（複数指定）
- `--group_labels` (default: `group1 group2`): 群名
- `--window` (default: 300): window 幅
- `--slide` (default: 200): slide 幅
- `--output` (default: `./matrix_out`): 出力ディレクトリ
- `--tmpdir`: 一時ディレクトリ


Output format (`<g1>_<g2>_<ctx>_matrix.tsv.gz`):
- 列順は `bedtools intersect -wa -wb` 由来の固定13列
- 実質内容: window座標 + site座標/群/sample/strand/meth/unmeth/coverage

Example (2 samples total by group):
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

Example (4 samples total by group):
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
- window 単位で GLMM を当てて p 値、delta、要約統計を推定。モデルを選択できるがbeta familyモデルsiteモードが最も正確性が高い。context, coverage, site数, replicate 数で段階フィルタ。必要に応じて `prefilter_delta` で高速事前除外も可能。


Options (major):
- `-i, --infile` (required): 入力 matrix
- `-o, --out_prefix` (required): 出力 prefix
- `-c, --context` (default: NULL): context フィルタ
- `--group1`, `--group2` (required): 比較群ラベル
- `--min_reps_g1`, `--min_reps_g2` (default: 2): 群ごとの最小 replicate
- `--family` (default: `beta`): `binom` または `beta`
- `--mode` (default: `site`): `aggregate` or `site`
- `--random_effect` (default: TRUE): `(1|sample)` を使用
- `--min_cov` (default: 0): site coverage フィルタ
- `--min_sites_win` (default: 0): window 内最小 site 数
- `--prefilter_delta` (default: 0): 事前デルタフィルタ
- `--workers` (default: 4): 並列 worker
- `--batches` (default: 50): 分割バッチ数
- `--max_globals_mb` (default: 1000): future global size
- `--seed` (default: 1): 乱数 seed

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
  --mode site \
  --min_reps_g1 2 --min_reps_g2 2 \
  --min_sites_win 3 \
  --min_cov 5 \
  --random_effect \
  --workers 8 --batches 200 --max_globals_mb 2000
```

## 3.5 `DMR_merge.R`

- window-level 有意シグナルを統合して DMR を構築。以下の５つのモードを実装している。

Supported merge modes:
- `Simes`: 隣接 window の p 値を Simes 法で統合し、複数の中程度シグナルをまとめて評価。
- `Stouffer`: 隣接 window の p 値を Z スコア化して統合し、広い領域の一貫した弱〜中シグナルを拾いやすい。
- `single_seed`: 強い単独 seed window を起点に extension して DMR を構築する、保守的で解釈しやすいモード。
- `multi_seed`: 有意 seed を複数含む領域を優先して連結するモードで、複数ピークを含む領域に強い。
- `hybrid_seed`: `single_seed` と `multi_seed` の考え方を併用し、感度と特異度のバランスを狙う実用的なデフォルト。

Options (major):
- `--windows` (required): GLMM window 結果ファイル（`*_fit_<family>_<mode>.tsv.gz`）
- `--out-prefix` (default: `results/dmr`): 出力 prefix
- `--merge-mode` (default: `hybrid_seed`): 上記5種のみ有効
- `--p-seed` (default: 0.05): seed 判定に使う p 閾値（全モードで使用）。
- `--p-extend` (default: 0.01): extension の許容 p 閾値（`Stouffer`, `single_seed`, `multi_seed`, `hybrid_seed` で使用。`Simes` では未使用）。
- `--max-gap-bp` (default: 200): 隣接 window を同一候補として連結する最大ギャップ（全モードで使用）。
- `--min-windows` (default: 2): DMR として採用する最小 window 数（全モードで使用）。
- `--merge-overlaps` / `--merge-overlaps-gap`: 検出後に同方向 DMR の重なり/近接を再マージ（`multi_seed` のみ有効）。
- `--post-filter`: 候補 DMR の品質フィルタを有効化（`Stouffer`, `single_seed`, `multi_seed`, `hybrid_seed` で有効、`Simes` では未使用）。
- `--min-median-p`, `--min-consistent-frac`: `--post-filter` 時の判定閾値（上記 post-filter 対応モードでのみ有効）。
- `--min-delta`: extension 時の最小効果量しきい値（`single_seed`, `multi_seed`, `hybrid_seed` で有効。`Simes`, `Stouffer` では未使用）。
- `--adaptive-delta`, `--adaptive-delta-method`, `--adaptive-delta-ratio`: delta しきい値をデータ依存で自動調整（`multi_seed` でのみ有効。`--adaptive-delta` を有効化した場合に使用）。
- `--max-p-degradation`, `--max-final-p`, `--min-strong-windows`: seed 拡張中および最終採択時の厳しさ制御（`single_seed`, `multi_seed` で有効）。
- `--trim-weak-edges`: DMR 両端の弱い window を削る（`single_seed`, `multi_seed` で有効）。
- `--min-dmr-length`: 最終 DMR 長の下限フィルタ（`single_seed`, `multi_seed` で有効）。
- `--max-median-p`: DMR 内 window の median p 上限（`single_seed`, `multi_seed` で有効）。
- `--seed-min-windows` (default: 1): multi-seed の seed 構成に必要な最小 window 数（`multi_seed` でのみ有効）。

Mode-wise option quick reference:

| Option group | Simes | Stouffer | single_seed | multi_seed | hybrid_seed |
|---|---|---|---|---|---|
| Base (`--p-seed`, `--max-gap-bp`, `--min-windows`) | yes | yes | yes | yes | yes |
| Extension threshold (`--p-extend`) | no | yes | yes | yes | yes |
| Post-filter (`--post-filter`, `--min-median-p`, `--min-consistent-frac`) | no | yes | yes | yes | yes |
| Delta threshold (`--min-delta`) | no | no | yes | yes | yes |
| Adaptive delta (`--adaptive-delta*`) | no | no | no | yes | no |
| Seed quality (`--max-p-degradation`, `--max-final-p`, `--min-strong-windows`) | no | no | yes | yes | no |
| Edge/length/median filters (`--trim-weak-edges`, `--min-dmr-length`, `--max-median-p`) | no | no | yes | yes | no |
| Overlap merge (`--merge-overlaps`, `--merge-overlaps-gap`) | no | no | no | yes | no |
| Multi-seed seed size (`--seed-min-windows`) | no | no | no | yes | no |

`--adaptive-delta*` = `--adaptive-delta`, `--adaptive-delta-method`, `--adaptive-delta-ratio`.

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

- site methylation （BinomTest 出力 TSV.gz）から bin ごとの平均 methylation を bigWig 化。

Options:
- `-i, --input` (required): binomtest 結果
- `-b, --binsize` (default: 50): bin 幅
- `-o, --output` (required): 出力 bigWig
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

Purpose:
- 複数 bigWig の bin 平均値から replicate 間分散を計算して bigWig 化。

Note:
- Python script。`python` で実行してください。

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
python make_binned_variance_bigwig.py \
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
