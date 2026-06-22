#!/usr/bin/env Rscript
# run_methylKit.R
suppressPackageStartupMessages({
  library(methylKit); library(data.table)
})

file.list <- list(
  "output_for_methylKit/sites_CG_forMethylKit_MT1.txt",
  "output_for_methylKit/sites_CG_forMethylKit_MT2.txt",
  "output_for_methylKit/sites_CG_forMethylKit_MT3.txt",
  "output_for_methylKit/sites_CG_forMethylKit_MT4.txt",
  "output_for_methylKit/sites_CG_forMethylKit_WT1.txt",
  "output_for_methylKit/sites_CG_forMethylKit_WT2.txt",
  "output_for_methylKit/sites_CG_forMethylKit_WT3.txt",
  "output_for_methylKit/sites_CG_forMethylKit_WT4.txt"
)
myobjDB <- methRead(
  file.list,
  sample.id = list("MT1","MT2","MT3","MT4",
                              "WT1","WT2","WT3","WT4"),
  assembly = "sim",
  treatment = c(1,1,1,1,0,0,0,0),
  context = "CpG",
  dbtype = "tabix",
  dbdir = "methylKit_DB"
)

filtered.myobj=filterByCoverage(myobjDB,lo.count=5,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

meth=unite(filtered.myobj, destrand=FALSE)

tiles <- tileMethylCounts(meth, win.size=300, step.size=200)
tiles_diff <- calculateDiffMeth(tiles, mc.cores=1)
sig_tiles <- getMethylDiff(tiles_diff, difference=0, qvalue=0.05)

write.table(getData(tiles_diff), file="output_for_methylKit/methylKit_diff.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(getData(sig_tiles), file="output_for_methylKit/methylKit_dmrs.tsv", sep="\t", quote=FALSE, row.names=FALSE)

sig_dt <- as.data.table(getData(sig_tiles))
if (nrow(sig_dt) > 0 && all(c("chr", "start", "end") %in% names(sig_dt))) {
  setorder(sig_dt, chr, start, end)
  merged_list <- list()
  cur <- sig_dt[1, .(chr, start, end)]
  if ("meth.diff" %in% names(sig_dt)) {
    cur$meth_diff_sum <- sig_dt$meth.diff[1]
    cur$meth_diff_n <- 1L
  }
  if (nrow(sig_dt) > 1) {
    for (i in 2:nrow(sig_dt)) {
      row <- sig_dt[i]
      if (row$chr == cur$chr && row$start <= (cur$end + 200)) {
        cur$end <- max(cur$end, row$end)
        if ("meth.diff" %in% names(sig_dt)) {
          cur$meth_diff_sum <- cur$meth_diff_sum + row$meth.diff
          cur$meth_diff_n <- cur$meth_diff_n + 1L
        }
      } else {
        if ("meth.diff" %in% names(sig_dt)) {
          cur$meth.diff <- cur$meth_diff_sum / cur$meth_diff_n
        }
        merged_list[[length(merged_list) + 1]] <- cur
        cur <- row[, .(chr, start, end)]
        if ("meth.diff" %in% names(sig_dt)) {
          cur$meth_diff_sum <- row$meth.diff
          cur$meth_diff_n <- 1L
        }
      }
    }
  }
  if ("meth.diff" %in% names(sig_dt)) {
    cur$meth.diff <- cur$meth_diff_sum / cur$meth_diff_n
  }
  merged_list[[length(merged_list) + 1]] <- cur
  merged_dmrs <- rbindlist(merged_list)
  if ("meth_diff_sum" %in% names(merged_dmrs)) {
    merged_dmrs[, `:=`(meth_diff_sum = NULL, meth_diff_n = NULL)]
  }
  fwrite(merged_dmrs, file="output_for_methylKit/methylKit_dmrs_merged.tsv", sep="\t")
}

