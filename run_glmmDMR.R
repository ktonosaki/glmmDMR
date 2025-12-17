#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(glmmTMB))
suppressMessages(library(parallel))
suppressMessages(library(furrr))
suppressMessages(library(data.table))

option_list <- list(
  make_option(c("-1", "--group1"), type="character", help="Group1 input directory"),
  make_option(c("-2", "--group2"), type="character", help="Group2 input directory"),
  make_option(c("-c", "--context"), type="character", default="CpG", help="Methylation context to analyze (CpG, CHG, CHH)"),
  make_option(c("-o", "--output_file"), type="character", help="Output TSV file (not ZIP file)"),
  make_option(c("-t", "--threads"), type="integer", default=8, help="Number of threads"),
  make_option("--parallel_engine", type="character", default="mclapply", help="mclapply or furrr"),
  make_option("--debug", action="store_true", default=FALSE, help="Enable debug output")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$group1) || is.null(opt$group2) || is.null(opt$output_file)) {
  cat("Missing required arguments\n")
  quit(status=1)
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

if (opt$debug) {
  cat("========== RUN CONFIGURATION ==========\n")
  cat("Context           :", opt$context, "\n")
  cat("Output file       :", opt$output_file, "\n")
  cat("Parallel engine   :", opt$parallel_engine, "\n")
  cat("Threads           :", opt$threads, "\n")
}

files_g1 <- get_file_list(opt$group1, opt$context)
files_g2 <- get_file_list(opt$group2, opt$context)
df1 <- load_data(files_g1, "group1")
df2 <- load_data(files_g2, "group2")
df <- bind_rows(df1, df2) %>%
  mutate(group = factor(group),
         sample = factor(sample),
         win_id = paste(chr, start, end, sep="_"))

fit_glmm <- function(sub) {
  if (nrow(sub) < 4) return(NULL)
  n_rep1 <- length(unique(sub$sample[sub$group == "group1"]))
  n_rep2 <- length(unique(sub$sample[sub$group == "group2"]))
  if (n_rep1 < 2 || n_rep2 < 2) return(NULL)
  
  sub <- make_offset(sub)
  sub$prop <- sub$meth / (sub$meth + sub$unmeth)
  sub <- sub[!is.nan(sub$prop), ]
  
  tryCatch({
      sub$adj_prop <- pmax(pmin(sub$prop, 0.99), 0.01)
      response <- "adj_prop"
      family <- glmmTMB::beta_family()
    
      formula_base <- paste(response, "~ group + (1 | sample)")
      formula_reduced <- paste(response, "~ 1 + (1 | sample)")
    
      full_model <- suppressWarnings(glmmTMB::glmmTMB(
        as.formula(formula_base),
        ziformula = ~0,
        data = sub,
        family = glmmTMB::beta_family()
      ))
      reduced_model <- suppressWarnings(glmmTMB::glmmTMB(
        as.formula(formula_reduced),
        ziformula = ~0,
        data = sub,
        family = glmmTMB::beta_family()
      ))

    lrt <- anova(reduced_model, full_model)
    pval <- tryCatch(lrt$`Pr(>Chisq)`[2], error=function(e) NA)
    aic_diff <- AIC(full_model) - AIC(reduced_model)
    bic_diff <- BIC(full_model) - BIC(reduced_model)
    
    delta <- diff(tapply(sub$prop, sub$group, mean, na.rm=TRUE))
    mean_rate1 <- mean(sub$prop[sub$group == "group1"], na.rm=TRUE)
    mean_rate2 <- mean(sub$prop[sub$group == "group2"], na.rm=TRUE)
    return(tibble(chr=sub$chr[1], start=sub$start[1], end=sub$end[1],
                  delta=delta, mean_rate1=mean_rate1, mean_rate2=mean_rate2,
                  pval=pval, aic_diff=aic_diff, bic_diff=bic_diff))
  }, error=function(e) NULL)
}

cat("Running GLMM with", opt$threads, "threads\n")
windows <- df %>% group_by(chr, start, end) %>% group_split()

if (opt$parallel_engine == "furrr") {
  plan(multisession, workers = opt$threads)
  results <- future_map_dfr(windows, fit_glmm, .progress=TRUE)
} else {
  results <- mclapply(windows, fit_glmm, mc.cores = opt$threads) %>% bind_rows()
}

results <- results %>% mutate(FDR = p.adjust(pval, method="BH"))
write_tsv(results, opt$output_file)

cat("Total windows to process:", length(windows), "\n")
cat("Output file will be written to:", opt$output_file, "\n")
cat("→ Written:", opt$output_file, "\n")
