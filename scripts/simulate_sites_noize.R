#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(GenomicRanges)
})

# ---------- options ----------
option_list <- list(
  make_option(c("-o","--out_dir"), type="character", default="results/site_window_sim",
              help="出力ディレクトリ"),
  make_option(c("--contexts"), type="character", default="CG,CHG,CHH",
              help="コンテキスト (カンマ区切り)"),
  make_option(c("--chr_len"), type="integer", default=5e6),
  make_option(c("--site_lambda_CG"), type="double", default=4e4),
  make_option(c("--site_lambda_CHG"), type="double", default=3e4),
  make_option(c("--site_lambda_CHH"), type="double", default=2e4),
  make_option(c("--block_frac"), type="double", default=0.10),
  make_option(c("--block_len_lo"), type="integer", default=500),
  make_option(c("--block_len_hi"), type="integer", default=5000),
  make_option(c("--delta_grid"), type="character", default="0.1,0.2,0.3"),
  make_option(c("--groups"), type="character", default="WT,MT"),
  make_option(c("--rep_per_group"), type="integer", default=3),
  make_option(c("--base_CG"), type="double", default=0.7),
  make_option(c("--base_CHG"), type="double", default=0.3),
  make_option(c("--base_CHH"), type="double", default=0.1),
  make_option(c("--dmr_site_prop"), type="double", default=1.0,
              help="各 DMR ブロックの中で本当にシフトするサイトの割合 (0〜1)"),
  make_option(c("--mu_site_sd"), type="double", default=0.0,
              help="サイトごとの平均メチル化率 mu に加えるランダム効果のSD"),
  make_option(c("--rho"), type="double", default=0.05, help="beta-binomialの過分散強さ（小さいほど強い）"),
  make_option(c("--zeta_CHH"), type="double", default=0.40, help="CHH用のゼロ過剰率"),
  make_option(c("--logcov_mu"), type="double", default=log(20)),
  make_option(c("--logcov_sd"), type="double", default=0.5),
  make_option(c("--miss_rate"), type="double", default=0.10, help="siteを欠測にする確率"),
  make_option(c("--min_cov"), type="integer", default=5, help="siteを有効とみなすカバレッジ閾値"),
  make_option(c("--seed"), type="integer", default=202509)
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(file.path(opt$out_dir, "tsv"), recursive=TRUE, showWarnings=FALSE)
set.seed(opt$seed)

# ---------- helpers ----------
rbeta_from_mean_conc <- function(mu, conc) {
  a <- mu*conc; b <- (1-mu)*conc
  stats::rbeta(1, max(a,1e-6), max(b,1e-6))
}
simulate_blocks <- function(L, frac, len_lo, len_hi, deltas) {
  n_try <- max(1, ceiling(frac * (L / ((len_lo+len_hi)/2))))
  st  <- pmax(1, sample.int(L, n_try))
  len <- sample(seq(len_lo, len_hi, by=10), n_try, replace=TRUE)
  en  <- pmin(L, st + len - 1L)
  d   <- sample(deltas, n_try, replace=TRUE) * sample(c(-1,1), n_try, replace=TRUE)
  GRanges("chr1", IRanges(st, en), dir=ifelse(d>0,"hyper","hypo"), delta=d)
}
simulate_sites <- function(L, lambda_perMb){
  n <- rpois(1, lambda_perMb * (L/1e6)); if (n==0) n <- 1
  pos <- sort(sample.int(L, n))
  data.table(chr="chr1", pos=pos)
}
gr_to_dt <- function(gr) {
  data.table(chr=as.character(seqnames(gr)),
             start=start(gr), end=end(gr),
             dir=as.character(mcols(gr)$dir),
             delta=as.numeric(mcols(gr)$delta))
}
label_truth_sites <- function(site_dt, blocks_gr) {
  if (length(blocks_gr)==0L || nrow(site_dt)==0L) {
    site_dt[, `:=`(truth=0L, dir=NA_character_)]
    return(site_dt)
  }
  sg  <- GRanges(site_dt$chr, IRanges(site_dt$pos, site_dt$pos))
  ov  <- findOverlaps(sg, blocks_gr, ignore.strand=TRUE)
  site_dt[, `:=`(truth=0L, dir=NA_character_)]
  if (length(ov)>0) {
    qh <- queryHits(ov); sh <- subjectHits(ov)
    site_dt$truth[qh] <- 1L
    site_dt$dir[qh]   <- as.character(mcols(blocks_gr)$dir[sh])
  }
  site_dt
}

# base rate / lambda per context
site_lambda <- list(CG=opt$site_lambda_CG, CHG=opt$site_lambda_CHG, CHH=opt$site_lambda_CHH)
base_rate   <- list(CG=opt$base_CG, CHG=opt$base_CHG, CHH=opt$base_CHH)
contexts    <- strsplit(opt$contexts, ",", fixed=TRUE)[[1]]
groups      <- strsplit(opt$groups, ",", fixed=TRUE)[[1]]

# DMR 周辺のノイズ抑制パラメータ
near_buffer  <- 200   # DMR から ±200 bp を「周辺」とみなす
near_factor  <- 0.3   # 周辺ではノイズを 0.3 倍にする

# ---------- simulate per context ----------
for (ctx in contexts) {
  message("=== Simulating sites: ", ctx, " ===")
  blocks <- simulate_blocks(opt$chr_len, opt$block_frac, opt$block_len_lo, opt$block_len_hi,
                            as.numeric(strsplit(opt$delta_grid,",",TRUE)[[1]]))
  sites  <- simulate_sites(opt$chr_len, lambda_perMb = site_lambda[[ctx]])
  # サンプルデザイン
  samples <- data.table(
    sample=sprintf("%s%02d", rep(groups, each=opt$rep_per_group),
                   sequence(rep(opt$rep_per_group, length(groups)))),
    group=rep(groups, each=opt$rep_per_group),
    replicate=rep(sequence(opt$rep_per_group), times=length(groups))
  )
  SS <- CJ(i=seq_len(nrow(sites)), s=seq_len(nrow(samples)))
  SS[, `:=`(chr=sites$chr[i], pos=sites$pos[i],
            sample=samples$sample[s], group=samples$group[s], replicate=samples$replicate[s])]
  SS[, context := ctx]
  
  # ブロック割り当て
  sg <- GRanges(SS$chr, IRanges(SS$pos, SS$pos))
  hit <- findOverlaps(sg, blocks, ignore.strand = TRUE)
  delta <- numeric(nrow(SS))
  in_block <- logical(nrow(SS))
  
  if (length(hit)>0) {
    tb <- as.data.table(hit)
    qh <- tb$queryHits; sh <- tb$subjectHits
    delta[qh]    <- mcols(blocks)$delta[sh]
    in_block[qh] <- TRUE
  }
  
  # --- DMR 内でも一部のサイトだけシフト（aggregate モデルを不利にする）---
  if (opt$dmr_site_prop < 1) {
    idx <- which(in_block)
    if (length(idx) > 0) {
      eff_flag <- rbinom(length(idx), size = 1, prob = opt$dmr_site_prop)
      delta[idx[eff_flag == 0]] <- 0  # このサイトは実質 non-DMR
    }
  }
  ## --- DMR 周辺サイトの判定（±near_buffer bp 拡張ブロックとのオーバーラップ）---
  if (length(blocks) > 0L) {
    blocks_buf <- blocks
    start(blocks_buf) <- pmax(1L, start(blocks_buf) - near_buffer)
    end(blocks_buf)   <- pmin(opt$chr_len, end(blocks_buf) + near_buffer)
    ov_near <- findOverlaps(sg, blocks_buf, ignore.strand = TRUE)
    is_near_dmr <- logical(nrow(SS))
    if (length(ov_near) > 0L) {
      is_near_dmr[queryHits(ov_near)] <- TRUE
    }
  } else {
    is_near_dmr <- logical(nrow(SS))
  }
  noise_factor <- ifelse(is_near_dmr, near_factor, 1.0)
  
  # 例：各 replicate に対してランダムオフセットを付加
  rep_sd <- 0.01  # 1% 程度のゆらぎ（調整可能）
  rep_effect <- rnorm(opt$rep_per_group * length(groups), mean = 0, sd = rep_sd)
  names(rep_effect) <- unique(SS$sample)
  
  # DMR 周辺では replicate ノイズを弱める
  mu <- base_rate[[ctx]] +
    delta * (SS$group==groups[2]) +
    rep_effect[SS$sample] * noise_factor
  
  # --- サイト間ヘテロを入れる（Beta モデル有利ゾーンを作る）---
  if (opt$mu_site_sd > 0) {
    mu <- mu + rnorm(length(mu), mean = 0, sd = opt$mu_site_sd)
  }
  
  mu <- pmin(pmax(mu, 1e-4), 1-1e-4)
  
  # coverage & 欠測（ここはそのまま）
  cov_rep_sd <- 0.1
  cov <- round(exp(rnorm(nrow(SS),
                         mean = opt$logcov_mu + rnorm(1, 0, cov_rep_sd),
                         sd   = opt$logcov_sd)))
  if (opt$miss_rate>0) cov[runif(length(cov)) < opt$miss_rate] <- 0L
  keep <- cov >= opt$min_cov
  SS <- SS[keep]; mu <- mu[keep]; cov <- cov[keep]; noise_factor <- noise_factor[keep]
  
  # 生成（CHHのみZI）
  use_zi <- identical(ctx,"CHH")
  meth <- unmeth <- rep(NA_integer_, nrow(SS))
  if (use_zi) {
    zi <- runif(nrow(SS)) < opt$zeta_CHH
    meth[zi]   <- 0L
    unmeth[zi] <- cov[zi]
  }
  
  # rho を site ごとに振り、DMR 周辺ではさらに弱める
  rho_site <- rlnorm(nrow(SS), meanlog=log(opt$rho), sdlog=0.15)
  rho_site <- rho_site * noise_factor          # 周辺では小さく（= conc 大きく = ばらつき↓）
  rho_site <- pmax(rho_site, 1e-6)            # 安全のため下限
  conc_site <- 1 / rho_site
  
  p_vec <- vapply(seq_along(mu),
                  function(i) rbeta_from_mean_conc(mu[i], conc_site[i]),
                  numeric(1))
  
  sel <- is.na(meth)
  if (any(sel)) {
    meth[sel]   <- rbinom(sum(sel), size=cov[sel], prob=p_vec[sel])
    unmeth[sel] <- pmax(cov[sel] - meth[sel], 0L)
  }
  
  DT <- data.table(chr=SS$chr, pos=SS$pos, sample=SS$sample, group=SS$group,
                   replicate=SS$replicate, context=ctx, meth=meth, unmeth=unmeth)
  # truth(site)
  DT_truth <- unique(DT[, .(chr,pos)])
  DT_truth <- label_truth_sites(DT_truth, blocks)
  setkey(DT, chr, pos, sample); setkey(DT_truth, chr, pos)
  DT <- DT_truth[DT]  # truth列をjoin
  setorder(DT, chr, pos, group, sample)
  
  fwrite(DT, file.path(opt$out_dir, "tsv", sprintf("sites_%s.tsv.gz", ctx)), sep="\t")
  fwrite(gr_to_dt(blocks), file.path(opt$out_dir, "tsv", sprintf("truth_blocks_%s.tsv.gz", ctx)), sep="\t")
  fwrite(unique(DT_truth[,.(chr,pos,truth,dir)]), file.path(opt$out_dir, "tsv", sprintf("truth_sites_%s.tsv.gz", ctx)), sep="\t")
}
message("Done: 01_simulate_sites.R")
