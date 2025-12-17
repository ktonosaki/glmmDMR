#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2) # Integrating ggplot2 for visualization
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

if (is.null(opt$windows)) stop("--windows is required")
if (opt$threads >= 0) data.table::setDTthreads(opt$threads)

#========================= Helper Functions ====================================
read_any <- function(path){
  message("[safe fread] ", path)
  tryCatch({
    if (grepl("\\.gz$", path)) {
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

calc_rho_adaptive <- function(D) {
  z <- z_from_p_twosided(D$p)
  rho <- rho_hat_lag1(z)
  return(rho)
}

stouffer_corr_adaptive <- function(pv, D) {
  rho <- calc_rho_adaptive(D)
  stouffer_corr(pv, rho=rho)
}

compare_methods <- function(sig_simes, sig_stouffer) {
  common <- nrow(fintersect(sig_simes, sig_stouffer))
  unique_simes <- nrow(fsetdiff(sig_simes, sig_stouffer))
  unique_stouffer <- nrow(fsetdiff(sig_stouffer, sig_simes))
  return(data.table(
    common = common,
    unique_simes = unique_simes,
    unique_stouffer = unique_stouffer,
    total_simes = nrow(sig_simes),
    total_stouffer = nrow(sig_stouffer)
  ))
}

plot_dmr_summary <- function(DMR, out_prefix) {
  summary_plot <- ggplot(DMR, aes(x=length_bp, y=delta_mean, color=direction)) +
    geom_point(alpha=0.7) +
    scale_x_log10() +
    labs(title="DMR Summary Plot",
         x="Region Length (bp, log scale)",
         y="Mean delta") +
    theme_minimal()
  ggsave(sprintf("%s.dmr_summary.pdf", out_prefix), plot=summary_plot, width=7, height=5)
}

# build_dmrs_dynamic \u95a2\u6570\u5185\u306e\u547c\u3073\u51fa\u3057\u3092\u66f4\u65b0
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
      new_p <- if (method == "simes") {
        p_simes(D$p[cand])
      } else {
        stouffer_corr_adaptive(D$p[cand], D)$p
      }
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

# DMR\u30d5\u30a1\u30a4\u30eb\u306e\u51fa\u529b\u51e6\u7406
DMR <- rbindlist(lapply(regions, pack_region))
plot_dmr_summary(DMR, opt$`out-prefix`)