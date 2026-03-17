#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input gzipped binomial test result (.tsv.gz)"),
  make_option(c("-b", "--binsize"), type="integer", default=50, help="Bin size in bp [default %default]"),
  make_option(c("-o", "--output"), type="character", help="Output bigWig file"),
  make_option("--genome", type="character", help="Genome .chrom.sizes file (chr<tab>size) for BigWig"),
  make_option("--context", type="character", default="CpG", help="Methylation context to use (CpG, CHG, CHH) [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate inputs
if (!file.exists(opt$input)) stop("Input file not found: ", opt$input)
if (!file.exists(opt$genome)) stop("Genome size file not found: ", opt$genome)
if (!opt$context %in% c("CpG", "CHG", "CHH")) stop("Invalid context: ", opt$context)

cat("[INFO] Reading input:", opt$input, "\n")
dt <- fread(cmd = paste("zcat", opt$input), header=TRUE)

# Filter by context
dt <- dt[context == opt$context]
if (nrow(dt) == 0) stop("No data found for context: ", opt$context)

# Calculate methylation rate
dt[, rate := meth / (meth + unmeth)]
dt <- dt[!is.na(rate) & is.finite(rate)]

# Bin assignment
dt[, bin := floor(pos / opt$binsize)]
dt[, bin_start := bin * opt$binsize]
dt[, bin_end := bin_start + opt$binsize]

# Summarize by bin
binned <- dt[, .(avg_rate = mean(rate, na.rm=TRUE)), by=.(chr, bin_start, bin_end)]
binned <- binned[!is.na(avg_rate)]

# Create GRanges
gr <- GRanges(
  seqnames = binned$chr,
  ranges = IRanges(start = binned$bin_start + 1, end = binned$bin_end),
  score = binned$avg_rate
)

# Load genome size
cat("[INFO] Loading genome size file:", opt$genome, "\n")
genome_sizes <- fread(opt$genome, header=FALSE, col.names=c("chr", "size"))
if (nrow(genome_sizes) == 0) stop("Empty genome size file")

valid_chr <- genome_sizes$chr

# Filter GRanges to chromosomes present in genome size file
cat("[INFO] Filtering chromosomes...\n")
orig_n <- length(gr)
gr <- gr[seqnames(gr) %in% valid_chr]
cat(sprintf("  Kept %d/%d ranges\n", length(gr), orig_n))

if (length(gr) == 0) stop("No ranges remaining after chromosome filtering")

# Set seqlevels and seqlengths with consistent ordering
seqlevels(gr) <- valid_chr
suppressWarnings(seqlengths(gr) <- setNames(genome_sizes$size, valid_chr))

# Trim ranges exceeding chromosome boundaries
gr <- trim(gr)

# Sort and export (suppress mitochondrial boundary warnings)
gr <- suppressWarnings(sort(gr))
cat("[INFO] Exporting BigWig to:", opt$output, "\n")
suppressWarnings(rtracklayer::export.bw(gr, opt$output))
cat("[INFO] Done.\n")
