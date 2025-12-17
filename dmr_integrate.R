#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

#========================= Options =============================================
option_list <- list(
  make_option(c("--windows"), type="character",
              help="GLMM window-level results (tsv/tsv.gz). Required."),
  make_option(c("--out-prefix"), type="character", default="results/dmr",
              help="Output prefix [default %default]"),
  make_option(c("--context"), type="character", default="CpG",
              help="Context label (for logs) [default %default]"),
  # p-value handling
  make_option(c("--p-side"), type="character", default="two-sided",
              help="Tail for p-value: two-sided/greater/less [default %default]"),
  # threshold merge parameters
  make_option(c("--p-seed"), type="numeric", default=0.05,
              help="Seed p threshold for starting a region [default %default]"),
  make_option(c("--p-extend"), type="numeric", default=0.05,
              help="Extend p threshold for adding windows [default %default]"),
  make_option(c("--max-gap-bp"), type="integer", default=200,
              help="Max gap in bp to continue merging [default %default]"),
  make_option(c("--max-bridge"), type="integer", default=1,
              help="Max number of non-significant bridges [default %default]"),
  make_option(c("--min-windows"), type="integer", default=2,
              help="Minimum #windows for a DMR [default %default]"),
  make_option(c("--min-length-bp"), type="integer", default=0,
              help="Filter DMR shorter than this bp [default %default]"),
  make_option(c("--min-delta"), type="numeric", default=0,
              help="Filter DMR with |delta_mean| < threshold [default %default]"),
  # dynamic merge parameters
  make_option(c("--merge-mode"), type="character", default="threshold",
              help="DMR merge strategy: threshold or dynamic [default %default]"),
  make_option(c("--p-dynamic"), type="numeric", default=0.05,
              help="Dynamic merging p threshold for combined p (Simes/Stouffer) [default %default]"),
  # FDR & threads
  make_option(c("--fdr"), type="numeric", default=0.05,
              help="FDR threshold (BH) [default %default]"),
  make_option(c("--threads"), type="integer", default=0,
              help="data.table threads (0=auto) [default %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

#opt$windows <- "/home/epigenome/dev_tools/glmm_dmrs/R_glmm_test/sim_from_empirical/sim_sites_fit_beta_aggregate.tsv.gz"
#opt$`p-seed`  <- 0.01
#opt$`p-extend` <- 0.05
#opt$`max-gap-bp`  <- 200
#opt$`max-bridge`  <- 0
#opt$`max-bridge`  <- 1
#opt$`min-windows`  <- 2
#opt$`--min-delta`  <- 0.1
#opt$`min-length-bp`  <- 300
#opt$`merge-mode`  <- "dynamic"
#opt$`merge-mode`  <- "threshold"



if (is.null(opt$windows)) stop("--windows is required")
if (opt$threads >= 0) data.table::setDTthreads(opt$threads)

#========================= Helper Functions ====================================
#read_any <- function(path){
#  if (grepl("\\.gz$", path)) fread(cmd = paste("zcat", shQuote(path))) else fread(path)
#}
#read_any <- function(path){
#  message("[safe fread] ", path)
#  tryCatch({
#    data.table::fread(path, fill = TRUE, showProgress = FALSE)
#  }, error = function(e){
#    message("[warning] fread failed, retrying with readLines() ...")
#    dat <- readLines(gzfile(path))
#    data.table::fread(text = dat, fill = TRUE)
#  })
#}

#read_any <- function(path){
#  message("[safe fread] ", path)
#  tryCatch({
#    if (grepl("\\.gz$", path)) {
#      # ---- 安全モード: zcat経由で読み込む（freadのzlibバグ回避） ----
#      data.table::fread(cmd = paste("zcat", shQuote(path)),
#                        fill = TRUE, showProgress = FALSE)
#    } else {
#      data.table::fread(path, fill = TRUE, showProgress = FALSE)
#    }
#  }, error = function(e){
#    message("[warning] fread failed, retrying with readLines() ...")
#    dat <- if (grepl("\\.gz$", path)) readLines(gzfile(path)) else readLines(path)
#    data.table::fread(text = dat, fill = TRUE)
#  })
#}

read_any <- function(path){
  message("[safe fread] ", path)
  tryCatch({
    if (grepl("\\.gz$", path)) {
      # ---- 完全安全モード: gzfile()でテキスト化してfreadに渡す ----
      con <- gzfile(path, "rt")
      on.exit(close(con))
      dat <- readLines(con)
      data.table::fread(text = dat, fill = TRUE, showProgress = FALSE)
    } else {
      data.table::fread(path, fill = TRUE, showProgress = FALSE)
    }
  }, error = function(e){
    message("[warning] fread(text) failed, retrying with base::read.table() ...")
    dat <- if (grepl("\\.gz$", path)) {
      read.table(gzfile(path), header = TRUE, sep = "\t", fill = TRUE, comment.char = "")
    } else {
      read.table(path, header = TRUE, sep = "\t", fill = TRUE, comment.char = "")
    }
    data.table::as.data.table(dat)
  })
}

tail_adjust <- function(p, delta, side){
  if (side == "two-sided") return(p)
  if (side == "greater")   return(ifelse(delta > 0, p/2, 1 - p/2))
  if (side == "less")      return(ifelse(delta < 0, p/2, 1 - p/2))
  p
}
p_simes <- function(pv){
  pv <- pv[is.finite(pv)]
  k <- length(pv); if (k == 0) return(1)
  pv <- sort(pv)
  min(1, min(pv * k / seq_len(k)))
}
z_from_p_twosided <- function(p){
  p2 <- pmin(1-1e-300, pmax(1e-300, p))
  qnorm(1 - p2/2)
}
rho_hat_lag1 <- function(z){
  n <- length(z); if (n <= 2) return(0)
  r <- suppressWarnings(cor(z[-n], z[-1], use="complete.obs"))
  r <- ifelse(is.finite(r), r, 0)
  pmin(0.99, pmax(0, r))  # clip to [0,0.99]
}
stouffer_corr <- function(pv, rho=0){
  m <- length(pv); if (m == 0) return(list(z=0, p=1))
  z <- z_from_p_twosided(pv)
  num <- sum(z)
  var_num <- (1 - rho)*m + rho*m^2
  z_c <- num / sqrt(var_num)
  p_c <- 2 * pnorm(-abs(z_c))
  list(z=z_c, p=p_c)
}

#========================= Load window data ====================================
cat("[load] windows: ", opt$windows, "\n")
win <- read_any(opt$windows)
#win <- fread(opt$windows)
req <- c("chr","start","end","p","delta")
if (!all(req %in% names(win))) stop("windows file must contain: ", paste(req, collapse=", "))

win[, p := tail_adjust(p, delta, opt$`p-side`)]
win <- win[is.finite(p) & p > 0 & p <= 1 & is.finite(delta)]
win[, direction := ifelse(delta >= 0, "hyper", "hypo")]
setorder(win, chr, start, end)

#========================= Merge: threshold ====================================
build_dmrs_threshold <- function(D){
  res <- list(); i <- 1L; n <- nrow(D)
  while (i <= n) {
    if (D$p[i] > opt$`p-seed`) { i <- i+1L; next }
    chr <- D$chr[i]; dir <- D$direction[i]
    run_idx <- i; include <- TRUE; bridge <- 0L; last_end <- D$end[i]; j <- i+1L
    while (j <= n && D$chr[j] == chr && D$direction[j] == dir) {
      gap <- D$start[j] - last_end
      if (gap > opt$`max-gap-bp`) break
      if (D$p[j] <= opt$`p-extend`) {
        run_idx <- c(run_idx, j); include <- c(include, TRUE); last_end <- max(last_end, D$end[j])
      } else if (bridge < opt$`max-bridge`) {
        run_idx <- c(run_idx, j); include <- c(include, FALSE); bridge <- bridge+1L; last_end <- max(last_end, D$end[j])
      } else break
      j <- j+1L
    }
    idx_keep <- run_idx[include]
    if (length(idx_keep) >= opt$`min-windows`) {
      region <- D[run_idx]
      res[[length(res)+1L]] <- list(chr=chr, start=min(region$start), end=max(region$end),
                                    n_windows=length(idx_keep), idx=idx_keep, direction=dir)
    }
    i <- j
  }
  res
}

#========================= Merge: dynamic (統合p<=閾値で継続) ==================
build_dmrs_dynamic <- function(D, method=c("simes","stouffer"), p_thresh=0.05){
  method <- match.arg(method)
  res <- list(); i <- 1L; n <- nrow(D)
  while (i <= n) {
    chr <- D$chr[i]; dir <- D$direction[i]
    if (D$p[i] > p_thresh) { i <- i+1L; next }
    run_idx <- i; j <- i + 1L
    while (j <= n && D$chr[j]==chr && D$direction[j]==dir) {
      gap <- D$start[j] - D$end[j-1]
      if (gap > opt$`max-gap-bp`) break
      cand <- c(run_idx, j)
      new_p <- if (method=="simes") p_simes(D$p[cand]) else stouffer_corr(D$p[cand], rho=0)$p
      if (new_p <= p_thresh) {
        run_idx <- cand; j <- j+1L
      } else break
    }
    if (length(run_idx) >= opt$`min-windows`) {
      res[[length(res)+1L]] <- list(chr=chr, start=min(D$start[run_idx]), end=max(D$end[run_idx]),
                                    n_windows=length(run_idx), idx=run_idx, direction=dir)
    }
    i <- max(run_idx) + 1L
  }
  res
}

#========================= Build regions =======================================
regions <- list()
if (opt$`merge-mode` == "threshold") {
  regions <- build_dmrs_threshold(win)
} else if (opt$`merge-mode` == "dynamic") {
  r1 <- build_dmrs_dynamic(win, "simes",    opt$`p-dynamic`)
  r2 <- build_dmrs_dynamic(win, "stouffer", opt$`p-dynamic`)
  regions <- c(r1, r2)
} else stop("Unknown --merge-mode")

if (length(regions) == 0L) {
  message("[done] No DMRs detected"); q(save="no", status=0)
}

#========================= Summarize each DMR ==================================
pack_region <- function(r){
  R <- win[r$idx]
  z <- z_from_p_twosided(R$p)
  rho <- rho_hat_lag1(z)
  st <- stouffer_corr(R$p, rho=rho)
  data.table(
    chr=r$chr,
    start=r$start,
    end=r$end,
    length_bp=r$end - r$start + 1L,
    n_windows=r$n_windows,
    direction=r$direction,
    delta_mean=mean(R$delta, na.rm=TRUE),
    p_simes=p_simes(R$p),
    z_stouffer=st$z,
    p_stouffer=st$p,
    rho_hat=rho
  )
}

DMR <- rbindlist(lapply(regions, pack_region))
if (opt$`min-length-bp` > 0) DMR <- DMR[length_bp >= opt$`min-length-bp`]
if (opt$`min-delta` > 0)     DMR <- DMR[abs(delta_mean) >= opt$`min-delta`]
DMR[, padj_BH_simes := p.adjust(p_simes, method="BH")]
DMR[, padj_BH_stouffer := p.adjust(p_stouffer, method="BH")]
DMR[, dmr_id := sprintf("DMR%06d", seq_len(.N))]

#========================= Outputs =============================================
dir.create(dirname(opt$`out-prefix`), showWarnings=FALSE, recursive=TRUE)
# すべてのDMRは常に保存（作業確認・再解析用）
fwrite(DMR, sprintf("%s.DMR.tsv.gz", opt$`out-prefix`), sep="\t")

if (opt$`merge-mode` == "dynamic") {
  #------ dynamic: Stouffer 有意のみを出力 ------------------------------------
  sig_stouffer <- DMR[padj_BH_stouffer < opt$fdr]
  
  if (nrow(sig_stouffer) > 0) {
    fwrite(sig_stouffer,
           sprintf("%s.DMR.sig.stouffer.tsv.gz", opt$`out-prefix`), sep="\t")
    
    # BED（score は -log10(p_stouffer)*100、name に p/FDR/Delta を入れる）
    for (dirn in c("hyper","hypo")) {
      B <- sig_stouffer[direction==dirn,
                        .(
#		          chr, start, end,
		          chr   = as.character(chr),
		          start = format(as.integer(start), scientific = FALSE, trim = TRUE),
			  end   = format(as.integer(end),   scientific = FALSE, trim = TRUE),
                          name = paste0(dmr_id,
                                        ",pZ:",   signif(p_stouffer,3),
                                        ",FDRZ:", signif(padj_BH_stouffer,3),
                                        ",Delta:",signif(delta_mean,3)),
                          score = round(-log10(p_stouffer)*100),
                          strand = ifelse(dirn=="hyper","+","-"))]
      if (nrow(B) > 0) {
        fwrite(B, sprintf("%s.%s.sig.bed", opt$`out-prefix`, dirn),
               sep="\t", col.names=FALSE, quote = FALSE)
      }
    }
  }
  
  # ログ要約（dynamic は Simes を “NA” 表示に）
  summary_log <- data.table(
    context      = opt$context,
    merge_mode   = opt$`merge-mode`,
    fdr          = opt$fdr,
    p_dynamic    = opt$`p-dynamic`,
    n_total      = nrow(DMR),
    n_sig_simes  = NA_integer_,
    n_sig_stouffer = nrow(sig_stouffer),
    n_sig_merged = NA_integer_,
    mean_length  = mean(DMR$length_bp),
    mean_delta   = mean(DMR$delta_mean),
    mean_rho_hat = mean(DMR$rho_hat)
  )
  fwrite(summary_log, sprintf("%s.summary.tsv", opt$`out-prefix`), sep="\t")
  
  message("[done/dynamic] DMRs=", nrow(DMR),
          " | sig(Stouffer)=", nrow(sig_stouffer))
  
} else {
  #------ threshold: 既存通りに Simes / Stouffer / merged を出力 --------------
  sig_simes     <- DMR[padj_BH_simes     < opt$fdr]
  sig_stouffer  <- DMR[padj_BH_stouffer  < opt$fdr]
  
  if (nrow(sig_simes) > 0)
    fwrite(sig_simes, sprintf("%s.DMR.sig.simes.tsv.gz", opt$`out-prefix`), sep="\t")
  if (nrow(sig_stouffer) > 0)
    fwrite(sig_stouffer, sprintf("%s.DMR.sig.stouffer.tsv.gz", opt$`out-prefix`), sep="\t")
  
  sig_all <- unique(rbindlist(list(sig_simes, sig_stouffer), fill=TRUE))
  if (nrow(sig_all) > 0) {
    fwrite(sig_all, sprintf("%s.DMR.merged.sig.tsv.gz", opt$`out-prefix`), sep="\t")
    for (dirn in c("hyper","hypo")) {
      B <- sig_all[direction==dirn,
                   .(
#                    chr, start, end,
                     chr   = as.character(chr),
                     start = format(as.integer(start), scientific = FALSE, trim = TRUE),
                     end   = format(as.integer(end),   scientific = FALSE, trim = TRUE),
                     name=paste0(dmr_id,
                                 ",pZ:",   signif(p_stouffer,3),
                                 ",FDRZ:", signif(padj_BH_stouffer,3),
                                 ",Delta:",signif(delta_mean,3)),
                     score=round(-log10(p_stouffer)*100),
                     strand=ifelse(dirn=="hyper","+","-"))]
      if (nrow(B)>0)
        fwrite(B, sprintf("%s.%s.sig.bed", opt$`out-prefix`, dirn),
               sep="\t", col.names=FALSE, quote = FALSE)
    }
  }
  
  summary_log <- data.table(
    context        = opt$context,
    merge_mode     = opt$`merge-mode`,
    fdr            = opt$fdr,
    p_seed         = opt$`p-seed`,
    p_extend       = opt$`p-extend`,
    n_total        = nrow(DMR),
    n_sig_simes    = nrow(sig_simes),
    n_sig_stouffer = nrow(sig_stouffer),
    n_sig_merged   = nrow(sig_all),
    mean_length    = mean(DMR$length_bp),
    mean_delta     = mean(DMR$delta_mean),
    mean_rho_hat   = mean(DMR$rho_hat)
  )
  fwrite(summary_log, sprintf("%s.summary.tsv", opt$`out-prefix`), sep="\t")
  
  message("[done/threshold] DMRs=", nrow(DMR),
          " | sig(Simes)=", nrow(sig_simes),
          " | sig(Stouffer)=", nrow(sig_stouffer),
          " | sig(merged)=", nrow(sig_all))
}
