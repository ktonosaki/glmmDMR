#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2) # ggplot2 integration for visualization
})

#========================= Options =============================================
option_list <- list(
  make_option(c("--windows"), type="character", help="GLMM window-level results (tsv/tsv.gz). Required."),
  make_option(c("--out-prefix"), type="character", default="results/dmr", help="Output prefix [default %default]"),
  make_option(c("--context"), type="character", default="CpG", help="Context label (for logs) [default %default]"),
  make_option(c("--merge-mode"), type="character", default="threshold", help="DMR merge strategy: threshold or dynamic [default %default]"),
  make_option(c("--p-seed"), type="numeric", default=0.05, help="Seed p threshold for starting a region [default %default]"),
  make_option(c("--p-extend"), type="numeric", default=0.05, help="Extend p threshold for adding windows [default %default]"),
  make_option(c("--max-gap-bp"), type="integer", default=200, help="Max gap in bp to continue merging [default %default]"),
  make_option(c("--max-bridge"), type="integer", default=1, help="Max number of non-significant bridges [default %default]"),
  make_option(c("--min-windows"), type="integer", default=2, help="Minimum #windows for a DMR [default %default]"),
  make_option(c("--min-length-bp"), type="integer", default=0, help="Filter DMR shorter than this bp [default %default]"),
  make_option(c("--p-dynamic"), type="numeric", default=0.05, help="Dynamic merging p threshold for combined p [default %default]"),
  make_option(c("--threads"), type="integer", default=0, help="data.table threads (0=auto) [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$windows)) stop("--windows is required")
if (opt$threads >= 0) data.table::setDTthreads(opt$threads)

#========================= Helper Functions ====================================
read_any <- function(path) {
  message("[DEBUG] Loading file: ", path)
  tryCatch({
    if (grepl("\\.gz$", path)) {
      con <- gzfile(path, "rt")
      on.exit(close(con))
      dat <- readLines(con)
      data.table::fread(text = dat)
    } else {
      data.table::fread(path)
    }
  }, error = function(e) {
    message("[ERROR] Failed to read file: ", path)
    stop(e)
  })
}

# Adjust tail for p-values
tail_adjust <- function(p, delta, side) {
  if (side == "two-sided") return(p)
  if (side == "greater") return(ifelse(delta > 0, p / 2, 1 - p / 2))
  if (side == "less") return(ifelse(delta < 0, p / 2, 1 - p / 2))
  return(p)
}

#========================= Load Data ===========================================
cat("[INFO] Loading window data\n")
win <- read_any(opt$windows)
required_cols <- c("chr", "start", "end", "p", "delta")
if (!all(required_cols %in% names(win))) {
  stop("Input file must contain columns: ", paste(required_cols, collapse = ", "))
}

cat("[INFO] Preprocessing windows\n")
win[, p := tail_adjust(p, delta, "two-sided")]
win <- win[is.finite(p) & p > 0 & p <= 1 & is.finite(delta)]
win[, direction := ifelse(delta >= 0, "hyper", "hypo")]
setorder(win, chr, start, end)

#========================= Merge Regions =======================================
build_dmrs_threshold <- function(D) {
  cat("[INFO] Using threshold-based merging\n")
  res <- list()
  i <- 1L; n <- nrow(D)
  while (i <= n) {
    if (D$p[i] > opt$`p-seed`) {
      i <- i + 1L
      next
    }
    chr <- D$chr[i]; dir <- D$direction[i]
    run_idx <- i
    last_end <- D$end[i]
    j <- i + 1L
    while (j <= n && D$chr[j] == chr && D$direction[j] == dir) {
      gap <- D$start[j] - last_end
      if (gap > opt$`max-gap-bp`) break
      if (D$p[j] <= opt$`p-extend`) {
        run_idx <- c(run_idx, j)
        last_end <- D$end[j]
      }
      j <- j + 1L
    }
    res[[length(res) + 1L]] <- list(chr = chr, start = min(D$start[run_idx]), end = max(D$end[run_idx]))
    i <- j
  }
  res
}

#========================= Execute Merge =======================================
if (opt$`merge-mode` == "threshold") {
  regions <- build_dmrs_threshold(win)
} else {
  stop("Unsupported merge mode: ", opt$`merge-mode`)
}

if (length(regions) == 0) {
  cat("[INFO] No regions detected\n")
  q(save = "no", status = 0)
}

cat("[INFO] Summarizing regions\n")
DMR <- data.table::rbindlist(lapply(regions, as.data.table))
DMR <- DMR[length_bp >= opt$`min-length-bp`]

#========================= Save Results ========================================
cat("[INFO] Saving results\n")
output_file <- sprintf("%s.DMR.tsv", opt$`out-prefix`)
data.table::fwrite(DMR, output_file, sep = "\t")