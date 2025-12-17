#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input gzipped binomial test result (.tsv.gz)"),
  make_option(c("-b", "--binsize"), type="integer", default=100, help="Bin size in bp [default %default]"),
  make_option(c("-o", "--output"), type="character", help="Output bigWig file"),
  make_option("--genome", type="character", help="Genome .chrom.sizes file (chr<tab>size) for BigWig"),
  make_option("--context", type="character", default="CpG", help="Methylation context to use (CpG, CHG, CHH) [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

#opt$input <- "/home/epigenome/rice_epigenome/parental_methylome/T2T_genome/05.binom_result/KxN_endosperm_rep1.genome1_binomtest_result.tsv.gz"
#opt$genome <- "/home/epigenome/reference/NIP-T2T/SNP_split/Kas_ref/Kas_n-masked_Mt_Pt.fa.txt"
#opt$binsize <- 50


cat("[INFO] Reading input:", opt$input, "\n")
dt <- fread(cmd = paste("zcat", opt$input), header=TRUE)
dt <- dt[context == opt$context]

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
gr
# Load genome size
cat("[INFO] Loading genome size file:", opt$genome, "\n")
genome_sizes <- fread(opt$genome, header=FALSE, col.names=c("chr", "size"))
# 明示的に名前を GRanges に合わせてセットする
chr_order <- as.character(seqlevels(gr))  # これが gr の順番と同じ
genome_sizes_matched <- genome_sizes[match(chr_order, genome_sizes$chr)]

# チェック
if (any(is.na(genome_sizes_matched$size))) {
  stop("Some chromosomes in GRanges not found in genome size file!")
}
# seqlengths に適用
seqlengths(gr) <- setNames(genome_sizes_matched$size, genome_sizes_matched$chr)

# Sort and export
gr <- sort(gr)
cat("[INFO] Exporting BigWig to:", opt$output, "\n")
rtracklayer::export.bw(gr, opt$output)
cat("[INFO] Done.\n")
