#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(data.table))
suppressMessages(library(furrr))
suppressMessages(library(future))

option_list <- list(
  make_option(c("-1", "--group1"), type = "character", help = "Group1 input directory or file"),
  make_option(c("-2", "--group2"), type = "character", help = "Group2 input directory or file"),
  make_option(c("-i", "--input_file"), type = "character", help = "Single input file with both groups (alternative to --group1/--group2)"),
  make_option(c("-o", "--output_file"), type = "character", help = "Output TSV file"),
  make_option(c("-c", "--context"), type = "character", default = "CG",
              help = "Methylation context to analyze (CG, CHG, CHH)"),
  make_option(c("--input_type"), type = "character", default = "sites",
              help = "Input data type: 'sites' (pos-based) or 'windows' (start-end based)"),
  make_option(c("--window_size"), type = "integer", default = 300,
              help = "Window size for site aggregation (only for --input_type sites)"),
  make_option(c("-t", "--threads"), type = "integer", default = 4,
              help = "Number of threads"),
  make_option(c("--parallel_engine"), type = "character", default = "mclapply", help = "mclapply or furrr"),
  make_option(c("--merge_dmrs"), action = "store_true", default = FALSE,
              help = "Merge significant windows into DMRs"),
  make_option(c("--fdr_threshold"), type = "numeric", default = 0.05,
              help = "FDR threshold for DMR detection"),
  make_option(c("--max_gap_bp"), type = "integer", default = 200,
              help = "Maximum gap (bp) between windows to merge into DMR"),
  make_option(c("--min_windows"), type = "integer", default = 2,
              help = "Minimum number of windows per DMR"),
  make_option(c("--debug"), action = "store_true", default = FALSE,
              help = "Enable debug output")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$output_file)) {
  cat("Missing required argument: --output_file\n")
  quit(status = 1)
}

if (is.null(opt$input_file) && (is.null(opt$group1) || is.null(opt$group2))) {
  cat("Must specify either --input_file OR both --group1 and --group2\n")
  quit(status = 1)
}

# Load site-level table
load_sites_data <- function(path) {
  if (opt$debug) cat("Loading sites data from:", path, "\n")
  read_tsv(path, show_col_types = FALSE) %>%
    filter(context == opt$context) %>%
    select(chr, pos, sample, group, context, meth, unmeth) %>%
    mutate(coverage = meth + unmeth)
}

# Load window-level table from one file
load_windows_data <- function(path) {
  if (opt$debug) cat("Loading windows data from:", path, "\n")
  read_tsv(path, show_col_types = FALSE,
           col_names = c("chr", "start", "end", "cytosine_chr", "cytosine_start", "cytosine_end",
                         "context", "group", "sample", "strand", "meth", "unmeth", "coverage")) %>%
    filter(context == opt$context) %>%
    select(chr, start, end, sample, group, context, meth, unmeth, coverage)
}

# Aggregate sites into fixed-size windows
aggregate_sites_to_windows <- function(sites, window_size) {
  if (opt$debug) cat("Aggregating sites to windows (size:", window_size, "bp)...\n")

  sites %>%
    mutate(
      win_start = floor(pos / window_size) * window_size,
      win_end = win_start + window_size
    ) %>%
    group_by(chr, win_start, win_end, sample, group, context) %>%
    summarise(
      meth = sum(meth, na.rm = TRUE),
      unmeth = sum(unmeth, na.rm = TRUE),
      coverage = sum(coverage, na.rm = TRUE),
      n_sites = n(),
      .groups = "drop"
    ) %>%
    rename(start = win_start, end = win_end) %>%
    mutate(win_id = paste(chr, start, end, sep = "_"))
}

get_file_list <- function(path, context) {
  list.files(path, pattern = paste0(context, ".*tsv.gz$"), full.names = TRUE)
}

load_data <- function(files, group_name) {
  bind_rows(lapply(files, function(f) {
    read_tsv(f, show_col_types = FALSE,
             col_names = c("chr", "start", "end", "cytosine_chr", "cytosine_start", "cytosine_end",
                           "context", "group", "sample", "strand", "meth", "unmeth", "coverage")) %>%
      mutate(group = group_name)
  }))
}

# Merge significant windows into DMR blocks
merge_windows_to_dmrs <- function(results, fdr_threshold, max_gap_bp, min_windows) {
  if (opt$debug) cat("Merging significant windows into DMRs...\n")

  sig_windows <- results %>%
    filter(FDR <= fdr_threshold, !is.na(pval)) %>%
    mutate(direction = ifelse(delta > 0, "hyper", "hypo")) %>%
    arrange(chr, start)

  if (nrow(sig_windows) == 0) {
    if (opt$debug) cat("No significant windows found.\n")
    return(tibble())
  }

  if (opt$debug) cat("Significant windows:", nrow(sig_windows), "\n")

  dmrs <- list()
  i <- 1
  n <- nrow(sig_windows)

  while (i <= n) {
    current_chr <- sig_windows$chr[i]
    current_dir <- sig_windows$direction[i]
    dmr_idx <- i

    j <- i + 1
    while (j <= n &&
           sig_windows$chr[j] == current_chr &&
           sig_windows$direction[j] == current_dir) {
      gap <- sig_windows$start[j] - sig_windows$end[j - 1]
      if (gap <= max_gap_bp) {
        dmr_idx <- c(dmr_idx, j)
        j <- j + 1
      } else {
        break
      }
    }

    if (length(dmr_idx) >= min_windows) {
      dmr_windows <- sig_windows[dmr_idx, ]

      combined_p <- tryCatch({
        chi_sq <- -2 * sum(log(dmr_windows$pval))
        df <- 2 * length(dmr_windows$pval)
        pchisq(chi_sq, df, lower.tail = FALSE)
      }, error = function(e) NA)

      dmrs[[length(dmrs) + 1]] <- tibble(
        chr = current_chr,
        start = min(dmr_windows$start),
        end = max(dmr_windows$end),
        n_windows = length(dmr_idx),
        direction = current_dir,
        delta_mean = mean(dmr_windows$delta, na.rm = TRUE),
        delta_median = median(dmr_windows$delta, na.rm = TRUE),
        mean_rate1 = mean(dmr_windows$mean_rate1, na.rm = TRUE),
        mean_rate2 = mean(dmr_windows$mean_rate2, na.rm = TRUE),
        combined_pval = combined_p,
        median_pval = median(dmr_windows$pval, na.rm = TRUE),
        min_pval = min(dmr_windows$pval, na.rm = TRUE)
      )
    }

    i <- max(dmr_idx) + 1
  }

  dmrs_df <- bind_rows(dmrs)

  if (nrow(dmrs_df) > 0) {
    dmrs_df <- dmrs_df %>%
      mutate(combined_FDR = p.adjust(combined_pval, method = "BH")) %>%
      arrange(chr, start)

    if (opt$debug) {
      cat("Total DMRs detected:", nrow(dmrs_df), "\n")
      cat("  Hypermethylated:", sum(dmrs_df$direction == "hyper"), "\n")
      cat("  Hypomethylated:", sum(dmrs_df$direction == "hypo"), "\n")
    }
  }

  return(dmrs_df)
}

if (opt$debug) cat("Reading input data...\n")

if (!is.null(opt$input_file)) {
  if (opt$input_type == "sites") {
    if (opt$debug) cat("Loading sites data and aggregating to windows...\n")
    sites <- load_sites_data(opt$input_file)
    df <- aggregate_sites_to_windows(sites, opt$window_size)
  } else {
    if (opt$debug) cat("Loading windows data...\n")
    df <- load_windows_data(opt$input_file)
  }
} else {
  if (opt$debug) cat("Loading data from directories...\n")
  files_g1 <- get_file_list(opt$group1, opt$context)
  files_g2 <- get_file_list(opt$group2, opt$context)
  df1 <- load_data(files_g1, "group1")
  df2 <- load_data(files_g2, "group2")
  df <- bind_rows(df1, df2)
}

# Common preprocessing
df <- df %>%
  mutate(
    group = factor(group),
    sample = factor(sample),
    win_id = paste(chr, start, end, sep = "_")
  )

unique_groups <- unique(df$group)
if (length(unique_groups) != 2) {
  cat("Error: Expected 2 groups, found", length(unique_groups), "\n")
  cat("Groups found:", paste(unique_groups, collapse = ", "), "\n")
  quit(status = 1)
}
GROUP1 <- as.character(unique_groups[1])
GROUP2 <- as.character(unique_groups[2])

if (opt$debug) cat("Groups detected:", GROUP1, "vs", GROUP2, "\n")

fit_fisher <- function(sub) {
  if (nrow(sub) < 2) {
    if (opt$debug) cat("Skipping window: not enough rows\n")
    return(NULL)
  }

  g_counts <- table(sub$group)
  if (length(g_counts) < 2 || any(g_counts < 1)) {
    if (opt$debug) cat("Skipping window: insufficient group data\n")
    return(NULL)
  }

  g1 <- sub %>% filter(group == GROUP1) %>%
    summarise(meth = sum(meth), unmeth = sum(unmeth))
  g2 <- sub %>% filter(group == GROUP2) %>%
    summarise(meth = sum(meth), unmeth = sum(unmeth))

  if (nrow(g1) == 0 || nrow(g2) == 0) return(NULL)

  mat <- matrix(c(g1$meth, g1$unmeth, g2$meth, g2$unmeth), nrow = 2)
  pval <- tryCatch(fisher.test(mat)$p.value, error = function(e) NA)

  sub$prop <- sub$meth / (sub$meth + sub$unmeth)
  mean_rate1 <- mean(sub$prop[sub$group == GROUP1], na.rm = TRUE)
  mean_rate2 <- mean(sub$prop[sub$group == GROUP2], na.rm = TRUE)
  delta <- mean_rate2 - mean_rate1

  tibble(
    chr = sub$chr[1], start = sub$start[1], end = sub$end[1],
    delta = delta, mean_rate1 = mean_rate1, mean_rate2 = mean_rate2,
    pval = pval
  )
}

windows <- df %>% group_by(chr, start, end) %>% group_split()

if (opt$parallel_engine == "furrr") {
  plan(multisession, workers = opt$threads)
  results <- future_map_dfr(windows, fit_fisher, .progress = TRUE)
} else {
  results_list <- parallel::mclapply(windows, fit_fisher, mc.cores = opt$threads)
  results <- dplyr::bind_rows(results_list[!sapply(results_list, is.null)])
}

results <- results %>% mutate(FDR = p.adjust(pval, method = "BH"))
write_tsv(results, opt$output_file)

if (opt$merge_dmrs) {
  dmrs <- merge_windows_to_dmrs(results, opt$fdr_threshold, opt$max_gap_bp, opt$min_windows)

  if (nrow(dmrs) > 0) {
    dmr_output <- sub("\\.tsv(\\.gz)?$", "_dmrs.tsv", opt$output_file)
    write_tsv(dmrs, dmr_output)
    if (opt$debug) cat("DMRs written to:", dmr_output, "\n")

    dmr_bed <- sub("\\.tsv(\\.gz)?$", "_dmrs.bed", opt$output_file)
    dmrs %>%
      select(chr, start, end, n_windows, direction) %>%
      write_tsv(dmr_bed, col_names = FALSE)
    if (opt$debug) cat("DMR BED file written to:", dmr_bed, "\n")
  } else {
    if (opt$debug) cat("No DMRs detected.\n")
  }
}

if (opt$debug) {
  cat("Total windows processed:", nrow(results), "\n")
  cat("Written to:", opt$output_file, "\n")
}
