# simulate_sites.R - CG Methylation Simulation

## 概要

`simulate_sites.R` は、GLMM-DMR用のテスト・ベンチマーク用に、現実的なCGメチル化データを**シミュレーション**するRスクリプトです。

### 主な機能

- **DMRブロック生成**：指定された効果量（delta）を持つDMR領域を自動生成
- **サイトレベルデータ生成**：各メチル化サイトについて、グループ間で異なるメチル化率を持つカウントデータを生成
- **複雑な分布モデル**：Beta-Binomial分布で過分散性を表現、zero-inflationにも対応
- **レプリケート間のノイズ**：現実的なサンプル変動とDMR周辺のノイズ低減を実装
- **グラウンドトゥルース出力**：検出精度評価用にDMRブロックと真のDMRサイトの情報を記録

---

## スクリプト構造

### 1. **ヘルパー関数**（Lines 39-82）

| 関数 | 役割 |
|-----|------|
| `rbeta_from_mean_conc(mu, conc)` | 平均と濃度パラメータからBeta分布のサンプリング |
| `simulate_blocks(L, frac, len_lo, len_hi, deltas)` | DMRブロック生成（GRanges形式） |
| `simulate_sites(L, lambda_perMb)` | サイト位置をPoisson分布でランダム生成 |
| `gr_to_dt(gr)` | GRanges → data.table変換 |
| `label_truth_sites(site_dt, blocks_gr)` | サイトがDMRに含まれるかフラグ付け |

### 2. **データ生成プロセス**（Lines 84-169）

```
DMRブロック生成
    ↓
サイト位置ランダムサンプリング
    ↓
DMRとサイトのオーバーラップ判定（重複時は絶対値最大のdelta採用）
    ↓
周辺領域判定（DMR±200bp）＆ノイズ調整
    ↓
グループ別メチル化率(mu)計算
    ↓
Beta-Binomial分布でメチル化カウント生成
    ↓
カバレッジフィルタ＆欠測処理
    ↓
出力ファイル生成
```

---

## CLIオプション詳細

### 基本設定

```bash
--out_dir DIR
  出力ディレクトリ（tsv/サブディレクトリに出力）
  [default: results/site_window_sim]

--chr_len N
  シミュレーション対象の染色体長（bp）
  [default: 5000000 (5Mbp)]

--seed N
  乱数シード（再現性のため）
  [default: 202509]
```

### サイト密度・メチル化率

```bash
--site_lambda LAMBDA
  サイト密度（1Mbあたりのサイト数）デフォルト40,000
  通常CG: 4e4, CHG: 3e4, CHH: 2e4
  [default: 4e4]

--base BASE_RATE
  ベース（野生型）メチル化率 (0-1)
  CG=0.7（70%メチル化）が標準
  [default: 0.7]
```

### DMRブロック設定

```bash
--block_frac FRAC
  シミュレーション領域のうちDMRで覆われる割合
  [default: 0.10] (10%)

--block_len_lo MIN_BP
  DMRブロックの最小長（bp）
  [default: 500]

--block_len_hi MAX_BP
  DMRブロックの最大長（bp）
  [default: 5000]

--delta_grid "0.1,0.2,0.3,..."
  DMRの効果量グリッド（メチル化率変化）
  例: 0.1 = MT群で±10%メチル化率が変わる
  複数値はカンマ区切り
  [default: "0.1,0.2,0.3"]
  
  注意：実際の効果量は符号がランダム（hyper/hypo）
```

### グループ・レプリケート設定

```bash
--groups "GROUP1,GROUP2"
  比較対象グループ（2群のみ固定）
  例: "WT,MT" や "Control,Treatment"
  [default: "WT,MT"]

--rep_per_group N
  各グループのレプリケート数
  [default: 3]
  
  例：--groups "WT,MT" --rep_per_group 4
  → WT01, WT02, WT03, WT04, MT01, MT02, MT03, MT04 生成
```

### DMR内パターン制御

```bash
--dmr_site_prop PROP
  各DMRブロック内で実際にメチル化率がシフトするサイトの割合
  [default: 1.0] (100% = すべてシフト)
  
  0.5に設定 → DMR内でも約50%のサイトのみ効果あり
  → aggregate法を不利にする現実的なシナリオ

--mu_site_sd SD
  サイト間ヘテロジェニティ（サイトごとのランダム効果）
  平均メチル化率muにN(0, SD)を加算
  [default: 0.0]
  
  0.05に設定 → ±5%程度のサイト固有のばらつき
```

### Beta-Binomial分布パラメータ

```bash
--rho RHO
  Beta-Binomial分布の過分散パラメータ（ρ）
  小さいほど過分散が強い（ばらつき大）
  [default: 0.05]
  
  0.05: 強い過分散（現実的）
  0.3: 過分散弱い（より規則的）
  
  内部演算: conc = 1 / rho
  rho が小さい → conc が大きい → Beta分布が尖る

--zeta_CHH RATE
  CHH文脈のみゼロ過剰率（現在CGは非ゼロ過剰固定）
  [default: 0.40] (40%)
  
  注意：現在scriptはCG固定なので本オプションは未使用
```

### カバレッジ・欠測制御

```bash
--logcov_mu MU
  カバレッジの対数平均
  [default: log(20) ≈ 2.996]
  
  exp(2.996) ≈ 20x coverage平均

--logcov_sd SD
  カバレッジの対数標準偏差
  [default: 0.5]
  
  レプリケート間でも異なる平均カバレッジを持つ

--miss_rate RATE
  各サイトが欠測（coverage=0）になる確率
  [default: 0.10] (10%)
  [推奨範囲: 0.05-0.20]

--min_cov MIN
  有効と判定する最小カバレッジ
  これ以下のサイトは除外
  [default: 5]
```

### 周辺領域ノイズ制御（内部固定）

```r
near_buffer  <- 200   # DMR±200bpを「周辺」判定
near_factor  <- 0.2   # 周辺ではノイズを0.2倍に低減
```

---

## 使用例

### 基本的な実行

```bash
Rscript simulate_sites.R
  # → デフォルトパラメータでCGデータを生成
  # → results/site_window_sim/tsv/ に出力
```

### ベンチマーク用（強いDMR効果）

```bash
Rscript simulate_sites.R \
  --out_dir benchmarks/strong_effect \
  --chr_len 10000000 \
  --site_lambda 20000 \
  --block_frac 0.05 \
  --delta_grid "0.3,0.5,0.7" \
  --block_len_lo 300 \
  --block_len_hi 2000 \
  --rep_per_group 4 \
  --rho 0.1 \
  --miss_rate 0.10
```

### 検査法の公平性テスト（弱いDMR）

```bash
Rscript simulate_sites.R \
  --out_dir benchmarks/weak_effect \
  --dmr_site_prop 0.5 \
  --delta_grid "0.05,0.1,0.15" \
  --rho 0.3 \
  --mu_site_sd 0.05 \
  --miss_rate 0.20
```

### 大規模シミュレーション

```bash
Rscript simulate_sites.R \
  --out_dir results/large_scale \
  --chr_len 50000000 \
  --site_lambda 30000 \
  --block_frac 0.08 \
  --rep_per_group 6 \
  --seed 123456
```

---

## 出力ファイル

### 1. `sites_CG.tsv.gz` - メチル化カウントデータ

**形式**: タブ区切り（gzip圧縮）

| 列名 | 型 | 説明 |
|------|-----|------|
| chr | char | 染色体 (固定: "chr1") |
| pos | integer | サイト位置（1-indexed） |
| sample | char | サンプルID (e.g., "WT01", "MT03") |
| group | char | グループ名 (e.g., "WT", "MT") |
| replicate | integer | レプリケート番号 (1, 2, 3, ...) |
| context | char | メチル化文脈 (固定: "CG") |
| meth | integer | メチル化リード数 |
| unmeth | integer | 非メチル化リード数 |
| truth | integer | DMR内フラグ (0=外部, 1=内部) |
| dir | char | DMR方向 (NA="非DMR", "hyper"=上昇型, "hypo"=低下型) |

**形式例**:
```
chr	pos	sample	group	replicate	context	meth	unmeth	truth	dir
chr1	1024	WT01	WT	1	CG	15	5	0	NA
chr1	2048	WT01	WT	1	CG	18	2	1	hyper
chr1	3072	MT01	MT	1	CG	8	12	1	hyper
```

### 2. `truth_blocks_CG.tsv.gz` - グラウンドトゥルース（DMRブロック）

**形式**: タブ区切り（gzip圧縮）

| 列名 | 型 | 説明 |
|------|-----|------|
| chr | char | 染色体 |
| start | integer | DMR開始位置（1-indexed） |
| end | integer | DMR終了位置（1-indexed, 含む） |
| dir | char | 効果方向 ("hyper"=WT→MT上昇, "hypo"=WT→MT低下) |
| delta | numeric | 効果量（メチル化率差） |

**形式例**:
```
chr	start	end	dir	delta
chr1	10000	15000	hyper	0.3
chr1	50000	52000	hypo	-0.1
chr1	100000	105000	hyper	0.2
```

### 3. `truth_sites_CG.tsv.gz` - グラウンドトゥルース（各サイト）

**形式**: タブ区切り（gzip圧縮）

| 列名 | 型 | 説明 |
|------|-----|------|
| chr | char | 染色体 |
| pos | integer | サイト位置 |
| truth | integer | DMR内フラグ (0 or 1) |
| dir | char | 属するDMRの方向（複数重複時は最大効果） |

**形式例**:
```
chr	pos	truth	dir
chr1	1024	0	NA
chr1	2048	1	hyper
chr1	12000	1	hyper
chr1	51500	1	hypo
```

---

## データ生成モデル

### メチル化率の階層モデル

```
mu_ij = base_rate + delta_i * I(group_j == "MT") + site_effect_i + rep_effect_j + noise_ij

where:
- base_rate = --base (WT群ベース)
- delta_i = サイト i が属するDMRの効果量（DMR外は0）
- I(.) = インジケータ関数
- site_effect_i = N(0, mu_site_sd) サイト間ヘテロ
- rep_effect_j = N(0, rep_sd * noise_factor) レプリケート効果
- noise_factor = 0.2 (DMR周辺), 1.0 (DMR外)
```

### メチル化カウント生成

```
p_{ij}k = rbeta(mu_ij, conc_ij)  ← Beta分布で確率抽出

meth_{ijl} ~ Binomial(cov_{ijl}, p_{ijk})

where:
- cov_{ijl} = coverage for site i, sample j, replicate l
- conc_ij = 1 / rho_ij  ← 濃度パラメータ
- rho_ij ~ LogNormal(log(--rho), 0.15) * noise_factor
```

### 複数DMR重複時の処理

サイト s が複数DMRと重複する場合：
```
delta_adopted = delta_{which.max( |delta| )}
                        ↑
                  絶対値が最大の効果量を採用
```

---

## パラメータ選択ガイド

### 実験条件によるチューニング

| シナリオ | 推奨設定 |
|--------|--------|
| **検査法の感度テスト** | `--delta_grid "0.05,0.1"` `--block_frac 0.02` 弱いDMRで性能評価 |
| **特異性（FDR）テスト** | `--delta_grid 0` または `--dmr_site_prop 0` ノイズのみで偽陽性率測定 |
| **現実的シナリオ** | `--block_frac 0.05-0.10` `--dmr_site_prop 0.5-0.8` データ品質低下を反映 |
| **理想的シナリオ** | `--block_frac 0.10` `--dmr_site_prop 1.0` `--miss_rate 0.05` 完全なDMR |
| **高ノイズ環境** | `--rho 0.3` `--miss_rate 0.30` `--mu_site_sd 0.1` 実験条件が悪い場合 |

### カバレッジ設計

```bash
# 低カバレッジ（RNA-Bisulfite-seqなど）
--logcov_mu 2.3 (≈ 10x) --logcov_sd 0.8 --miss_rate 0.20

# 標準カバレッジ（Whole-genome Bisulfiteなど）
--logcov_mu 2.996 (≈ 20x) --logcov_sd 0.5 --miss_rate 0.10

# 高カバレッジ（ターゲット化Bisulfiteなど）
--logcov_mu 3.9 (≈ 50x) --logcov_sd 0.3 --miss_rate 0.05
```

---

## 実行上の注意

### メモリ使用量

- `chr_len=5M + site_lambda=4e4` → 約200サイト × 8レプリケート = 1600行データテーブル
- `chr_len=50M + site_lambda=3e4` → 約1500サイト × 8レプリケート = 12000行
- 各セルが ~200byte → 2.4MB（メモリ効率的）

### 実行時間

典型的な環境（R 4.0+、Linux環境）：
- デフォルト設定（5Mbp）: 1-3秒
- 大規模（50Mbp）: 10-30秒
- 非常に大規模（500Mbp）: 数分

### 乱数再現性

同じ `--seed` を使用すれば、異なる実行でも同一データが再現されます。
ベンチマーク比較時に有用です。

---

## データの検証メソッド

### 生成データの基本統計

```r
library(data.table)
DT <- fread("results/site_window_sim/tsv/sites_CG.tsv.gz")

# 全体メチル化率
DT[, .(meth_rate = sum(meth) / (sum(meth) + sum(unmeth))),
   by = .(group)]

# グループ別メチル化率
DT[truth == 1, .(mean_meth = mean(meth / (meth + unmeth))),
   by = .(group, dir)]

# カバレッジ分布
DT[, cov := meth + unmeth]
DT[, .(mean_cov = mean(cov), sd_cov = sd(cov)), by = group]
```

### グラウンドトゥルース検証

```r
truth_sites <- fread("results/site_window_sim/tsv/truth_sites_CG.tsv.gz")
truth_blocks <- fread("results/site_window_sim/tsv/truth_blocks_CG.tsv.gz")

# DMR内のサイト数確認
truth_sites[truth == 1, .N]

# DMRブロックの統計
truth_blocks[, .(n_blocks = .N, mean_len = mean(end - start + 1),
                 n_hyper = sum(dir == "hyper"), n_hypo = sum(dir == "hypo"))]
```

---

## トラブルシューティング

### 空のDMRが生成されない場合

```bash
# --block_fracを増加させる
--block_frac 0.20  # 20%に増加

# または短いブロック長を設定
--block_len_lo 100 --block_len_hi 500
```

### 出力ファイルがない場合

1. 出力ディレクトリへの書き込み権限確認
2. ディスク容量確認
3. Rのエラーメッセージを確認

```bash
Rscript simulate_sites.R --out_dir /tmp/test 2>&1
```

### データが期待と異なる場合

1. `--seed`を除外して報告（再現性依存）
2. パラメータ設定を明記して確認
3. 出力の最初と最後の数行をチェック

---

## 参考・応用

このシミュレーションは以下の目的に使用できます：

- ✅ GLMM-DMRの感度・特異性評価
- ✅ 異なる統計手法の比較
- ✅ サンプルサイズ設計
- ✅ 多重比較補正の効果検証
- ✅ 検出力分析

---

## ライセンス・引用

本スクリプトはglmmDMRプロジェクトの一部です。
使用時は適切に引用してください。
