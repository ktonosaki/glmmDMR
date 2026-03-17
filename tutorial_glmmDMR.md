# glmmDMR Tutorial

このページは、glmmDMRパイプラインを実データで動かすときの詳細手順と、主要オプションの意味をまとめた実行ガイドです。
READMEは概要、こちらは実行デモとパラメータ設計の解説を担当します。

## 1. Pipeline overview

標準フローは次の順です。

1. Read preprocessing and methylation extraction (upstream)
2. summarize_extractor.py で site 集計
3. BinomTest.py で site ごとの有意性フィルタ
4. make_binned_methylation_bigwig.R と make_binned_variance_bigwig.R で可視化トラック生成 (任意)
5. prepare_matrix.sh で 2群比較用 window matrix を作成
6. run_glmmDMR.R で window ごとの GLMM を実行
7. DMR_merge.R で window を DMR に統合
8. 必要に応じて他手法比較や genomic context annotation

## 2. Prerequisites

必須ソフト (本リポジトリの主要スクリプト実行):
- Python packages: pandas, numpy, scipy, statsmodels, pyBigWig, tqdm
- R packages: optparse, data.table, glmmTMB, future, future.apply, GenomicRanges, rtracklayer
- CLI tools: bedtools, samtools, zcat

上流工程でよく使うツール:
- fastp
- SRA Toolkit (fasterq-dump)
- Bismark + bowtie2
- deepTools (bigwigAverage など)

推奨環境変数:
- TMPDIR を大容量ストレージに設定

例:
export TMPDIR=/path/to/large_storage/tmp
mkdir -p "$TMPDIR"

## 3. Step-by-step demo

### Step 1. Read preprocessing and methylation extraction

目的:
- 生 read を QC/trim し、Bismark extractor 出力を作る

最小例 (paired-end):
fasterq-dump SRRXXXXXXX --split-files --threads 8 --outdir raw
pigz -p 8 raw/SRRXXXXXXX_1.fastq
pigz -p 8 raw/SRRXXXXXXX_2.fastq

fastp \
  -i raw/SRRXXXXXXX_1.fastq.gz \
  -I raw/SRRXXXXXXX_2.fastq.gz \
  -o clean/sample_1.trimmed.fq.gz \
  -O clean/sample_2.trimmed.fq.gz \
  --thread 8

bismark --bowtie2 -p 8 -o align /path/to/bismark_genome \
  -1 clean/sample_1.trimmed.fq.gz \
  -2 clean/sample_2.trimmed.fq.gz

deduplicate_bismark --paired --output_dir align --bam align/sample_bismark_bt2_pe.bam

bismark_methylation_extractor --paired-end --gzip --parallel 8 \
  --genome_folder /path/to/bismark_genome \
  --output calls/sample \
  align/sample_bismark_bt2_pe.deduplicated.bam

成果物:
- calls/sample 以下に context 別 extractor ファイル

### Step 2. summarize_extractor.py

目的:
- extractor 出力を site 単位の集計テーブルに変換

実行例:
python summarize_extractor.py \
  -i calls/sample \
  -o matrix/sample_summarized_output.tsv.gz \
  --threads 4

主要オプション:
- -i, --input: extractor 出力ディレクトリ
- -o, --output: 出力 TSV.GZ
- --threads: 並列数
- --keep-temp: 中間ファイルを保持

出力列:
- chr, pos, strand, meth, unmeth, context

### Step 3. BinomTest.py

目的:
- site 単位で binomial test と FDR 補正を行い、低信頼シグナルを抑制

実行例:
python BinomTest.py \
  -i matrix/sample_summarized_output.tsv.gz \
  -o binom/sample_binomtest_result.tsv.gz \
  --nonconv_chr chloroplast \
  --min_coverage 5 \
  --fdr_threshold 0.05 \
  --threads 4

主要オプション:
- --null_prob: 帰無確率を固定指定 (未指定時は nonconv_chr 推定 or 0.5)
- --fdr_threshold: 有意判定の閾値
- --min_coverage: 最低カバレッジ
- --nonconv_chr: 変換効率推定に使う染色体名
- --threads: 並列数

動作メモ:
- 非有意 site は meth を 0 に置換して出力 (weighted-level downstream 向け)

### Step 4. bigWig generation (optional)

4-1. 平均メチル化 bigWig

実行例:
Rscript make_binned_methylation_bigwig.R \
  -i binom/sample_binomtest_result.tsv.gz \
  -b 50 \
  -o wig/sample_CpG.bw \
  --genome /path/to/TAIR10.chromInfo \
  --context CpG

主要オプション:
- -i: BinomTest 出力
- -b, --binsize: bin 幅
- -o: bigWig 出力
- --genome: chrom.sizes
- --context: CpG, CHG, CHH

4-2. replicate 分散 bigWig

実行例:
python make_binned_variance_bigwig.R \
  --inputs wig/rep1_CpG.bw wig/rep2_CpG.bw wig/rep3_CpG.bw wig/rep4_CpG.bw \
  --output wig/group_CpG.variance.bw \
  --bin-size 200 \
  --min-tracks 2 \
  --norm none

主要オプション:
- --inputs: 対象 bigWig 一覧
- --output: 分散 bigWig
- --bin-size: 集約 bin
- --min-tracks: 分散計算に必要な最小トラック数
- --norm: none または log2p1

### Step 5. prepare_matrix.sh

目的:
- 2群比較用の sliding-window matrix 作成

実行例:
bash prepare_matrix.sh \
  --fasta /path/to/TAIR10.fasta \
  --group1 binom/WT_1_binomtest_result.tsv.gz binom/WT_2_binomtest_result.tsv.gz \
  --group2 binom/MT_1_binomtest_result.tsv.gz binom/MT_2_binomtest_result.tsv.gz \
  --group_labels WT MT \
  --window 500 \
  --slide 300 \
  --output prep_out

主要オプション:
- --fasta: FASTA または FAI
- --group1, --group2: 比較2群の site テーブル
- --group_labels: 出力に使う群名
- --window: window 幅
- --slide: スライド幅
- --output: 出力ディレクトリ
- --tmpdir: 一時領域

注意:
- window を小さくすると解像度は上がるが計算量が増加
- slide を小さくすると連続性は上がるが相関が強くなる

### Step 6. run_glmmDMR.R

目的:
- window ごとに GLMM を当てて p と delta を推定

実行例:
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

主要オプション:
- -i, -o: 入力 matrix と出力 prefix
- --group1, --group2: 群ラベル
- --family: binom or beta
- --mode: aggregate or site
- --min_reps_g1, --min_reps_g2: 群ごとの最低 replicate 数
- --min_cov: site 最低 coverage
- --min_sites_win: window 内最低 site 数
- --prefilter_delta: 低効果量 window の事前除外
- --random_effect: sample ランダム効果の有無
- --workers, --batches, --max_globals_mb: 並列実行制御

出力:
- *_fit_family_mode.tsv.gz (p, delta, mean rates, AIC/BIC差など)

### Step 7. DMR_merge.R

目的:
- window-level signal を連結し DMR として統合

実行例:
Rscript DMR_merge.R \
  --windows glmm_out/WT_MT_CpG_fit_beta_aggregate.tsv.gz \
  --out-prefix dmr_out/WT_MT_CpG \
  --merge-mode stouffer_multi_seed \
  --p-seed 0.05 \
  --p-extend 0.05 \
  --min-windows 1

主要オプション (共通):
- --windows: GLMM出力
- --out-prefix: 出力接頭辞
- --merge-mode: 統合法
- --p-seed: シード判定閾値
- --p-extend: 拡張判定閾値
- --max-gap-bp: 窓間最大ギャップ
- --min-windows: DMR 最低窓数

主要オプション (品質制御):
- --min-delta, --adaptive-delta 系
- --max-p-degradation
- --max-final-p
- --min-strong-windows
- --trim-weak-edges
- --merge-overlaps, --merge-overlaps-gap

主な merge-mode:
- single_seed
- stouffer_multi_seed
- hybrid_seed
- combined
- simes_extended
- fdr_filter
- hierarchical
- weighted

### Step 8. Optional downstream analyses

例:
- DSS, methylKit, DMRfinder との比較
- promoter / gene body / TE / intergenic との overlap 集計
- context 別の DMR 分布比較

## 4. Recommended parameter presets

少数 replicate (n=2 vs n=2) の例:
- run_glmmDMR.R: --family beta --mode aggregate --min_reps_g1 2 --min_reps_g2 2 --min_sites_win 1 --min_cov 5
- DMR_merge.R: --merge-mode stouffer_multi_seed --p-seed 0.05 --p-extend 0.05 --min-windows 1

中規模 replicate (n=4 vs n=4) の例:
- run_glmmDMR.R: --family beta --mode aggregate --min_reps_g1 4 --min_reps_g2 4
- DMR_merge.R: --merge-mode stouffer_multi_seed または hybrid_seed

## 5. Troubleshooting

よくある問題:

1) No windows after filtering
- 原因: --min_cov, --min_sites_win, replicate 条件が厳しすぎる
- 対応: 閾値を段階的に緩める

2) GLMM が遅い / メモリ不足
- 対応: --workers を下げる、--batches を増やす、TMPDIR を高速かつ容量十分な場所へ

3) bigWig 生成で chromosome mismatch
- 原因: 入力 chr 名と chromInfo の不一致
- 対応: chr naming を統一してから再実行

4) DMR が細切れになる
- 対応: --max-gap-bp を拡大、merge-mode を再検討、p-extend と min-delta を調整

## 6. Reproducibility checklist

- 参照ゲノムと chromInfo のバージョンを固定
- 乱数 seed を固定 (run_glmmDMR.R の --seed)
- 依存パッケージのバージョンを記録
- 実行コマンドをログ保存
