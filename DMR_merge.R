#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(stats)
})

#========================= Options =============================================
option_list <- list(
  make_option(c("--windows"), type="character",
              help="GLMM window-level results (tsv/tsv.gz). Required."),
  make_option(c("--out-prefix"), type="character", default="results/dmr",
              help="Output prefix for results [default %default]"),
  make_option(c("--merge-mode"), type="character", default="hybrid_seed",
              help="DMR merge strategy: Simes, Stouffer, single_seed, multi_seed, hybrid_seed [default %default]"),
  make_option(c("--p-seed"), type="numeric", default=0.05,
              help="Seed p-value threshold for starting a DMR [default %default]"),
  make_option(c("--p-extend"), type="numeric", default=0.01,
              help="Extend p-value threshold for growing a DMR [default %default]"),
  make_option(c("--max-gap-bp"), type="integer", default=200,
              help="Max gap in bp to continue merging [default %default]"),
  make_option(c("--merge-overlaps"), action="store_true", default=FALSE,
              help="Merge overlapping DMRs after detection [default %default]"),
  make_option(c("--merge-overlaps-gap"), type="integer", default=0,
              help="Max gap (bp) to merge overlapping DMRs [default %default]"),
  make_option(c("--min-windows"), type="integer", default=2,
              help="Minimum windows for a DMR [default %default]"),
  make_option(c("--post-filter"), action="store_true", default=FALSE,
              help="Enable post-filtering to reduce false positives [default %default]"),
  make_option(c("--min-median-p"), type="numeric", default=0.01,
              help="Minimum median p-value for DMRs (post-filter) [default %default]"),
  make_option(c("--min-consistent-frac"), type="numeric", default=0.5,
              help="Minimum fraction of windows with p < p-seed (post-filter) [default %default]"),
  make_option(c("--min-delta"), type="numeric", default=0,
              help="Minimum absolute delta for extension (effect size filter) [default %default]"),
  make_option(c("--max-p-degradation"), type="numeric", default=1.2,
              help="Maximum allowed p-value degradation during extension (1.0=no degradation) [default %default]"),
  make_option(c("--max-final-p"), type="numeric", default=1.0,
              help="Maximum final combined p-value for DMRs (strict filter) [default %default]"),
  make_option(c("--min-strong-windows"), type="numeric", default=0.5,
              help="Minimum fraction of windows with p < p-seed in final DMR [default %default]"),
  make_option(c("--trim-weak-edges"), action="store_true", default=FALSE,
              help="Trim weak windows (p > p-extend) from DMR edges [default %default]"),
  make_option(c("--min-dmr-length"), type="integer", default=0,
              help="Minimum DMR length in bp (0=no filter) [default %default]"),
  make_option(c("--max-median-p"), type="numeric", default=1.0,
              help="Maximum median p-value for final DMRs [default %default]"),
  make_option(c("--seed-min-windows"), type="integer", default=1,
              help="Minimum windows for multi-seed Stouffer method (1=single+multi, 2=multi only) [default %default]"),
  make_option(c("--adaptive-delta"), action="store_true", default=FALSE,
              help="Use adaptive delta threshold based on candidate DMR distribution [default %default]"),
  make_option(c("--adaptive-delta-method"), type="character", default="median_ratio",
              help="Adaptive delta method: median_ratio, q50, q25, q10, mad [default %default]"),
  make_option(c("--adaptive-delta-ratio"), type="numeric", default=0.6,
              help="Ratio for median_ratio method (0.6 = 60%% of median) [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$windows)) stop("--windows argument is required.")

#========================= Helper Functions ====================================
load_data <- function(path) {
  message("[INFO] Loading data: ", path)
  dt <- fread(path)
  stopifnot(c("chr", "start", "end", "p", "delta") %in% names(dt))
  dt <- dt[order(chr, start)]
  dt[, direction := ifelse(delta >= 0, "hyper", "hypo")]
  return(dt[is.finite(p) & p > 0 & p <= 1 & is.finite(delta)])
}

# BEDファイル出力
write_bed <- function(dmrs, out_path) {
  bed_data <- dmrs[, .(
    chr,
    start = as.integer(start - 1L),
    end = as.integer(end),
    score = as.integer(n_windows),
    direction
  )]
  fwrite(bed_data, out_path, sep = "\t", col.names = FALSE)
}

p_simes <- function(p_values) {
  p_values <- p_values[is.finite(p_values)]
  k <- length(p_values)
  if (k == 0) return(1)
  p_sorted <- sort(p_values)
  min(1, min(p_sorted * k / seq_len(k)))
}

p_stouffer <- function(p_values) {
  p_values <- p_values[is.finite(p_values)]
  k <- length(p_values)
  if (k == 0) return(1)
  z_scores <- qnorm(1 - p_values)
  combined_z <- sum(z_scores) / sqrt(k)
  2 * pnorm(-abs(combined_z))
}

# DMRごとのdelta計算
add_delta_summary <- function(dmrs, win) {
  if (nrow(dmrs) == 0) return(dmrs)
  dmrs_with_delta <- lapply(1:nrow(dmrs), function(i) {
    row <- dmrs[i, ]
    overlapping_windows <- win[chr == row$chr &
                                  start >= row$start &
                                  end <= row$end]
    delta_mean <- mean(overlapping_windows$delta, na.rm = TRUE)
    delta_max <- max(overlapping_windows$delta, na.rm = TRUE)
    row[, `:=`(delta_mean = delta_mean, delta_max = delta_max)]
    return(row)
  })
  return(data.table::rbindlist(dmrs_with_delta))
}

# 適応的delta閾値の計算: 候補DMRのdelta分布から閾値を決定
calculate_adaptive_delta_threshold <- function(dmrs, win, method = "median_ratio", ratio = 0.6) {
  if (nrow(dmrs) == 0) {
    message("[WARN] No DMRs found for adaptive threshold calculation. Using default min-delta.")
    return(0)
  }
  
  # 各DMRに含まれるウィンドウのdelta値を収集
  all_deltas <- unlist(lapply(1:nrow(dmrs), function(i) {
    row <- dmrs[i, ]
    overlapping_windows <- win[chr == row$chr &
                                  start >= row$start &
                                  end <= row$end]
    abs(overlapping_windows$delta)
  }))
  
  all_deltas <- all_deltas[is.finite(all_deltas) & all_deltas > 0]
  
  if (length(all_deltas) == 0) {
    message("[WARN] No valid delta values found. Using default min-delta.")
    return(0)
  }
  
  # 統計量の計算
  median_delta <- median(all_deltas, na.rm = TRUE)
  q25_delta <- quantile(all_deltas, 0.25, na.rm = TRUE)
  q10_delta <- quantile(all_deltas, 0.10, na.rm = TRUE)
  mad_delta <- mad(all_deltas, na.rm = TRUE)  # Median Absolute Deviation
  
  # 閾値決定方法
  threshold <- switch(method,
    "median_ratio" = median_delta * ratio,     # 中央値の60%など
    "q50" = median_delta,                      # 中央値（50パーセンタイル）
    "q25" = q25_delta,                         # 第1四分位点
    "q10" = q10_delta,                         # 10パーセンタイル
    "mad" = median_delta - mad_delta,          # 中央値 - MAD
    median_delta * ratio                       # デフォルト
  )
  
  message(sprintf("[INFO] Adaptive delta threshold: %.4f (method=%s, median=%.4f, Q25=%.4f, Q10=%.4f)",
                  threshold, method, median_delta, q25_delta, q10_delta))
  
  return(max(0, threshold))  # 負の値を防ぐ
}

# DMR境界の最適化: 両端の弱いウィンドウを除去
trim_dmr_edges <- function(dmrs, win, p_extend) {
  if (nrow(dmrs) == 0) return(dmrs)
  
  message("[INFO] Trimming weak edges (p > ", p_extend, ")")
  
  trimmed <- lapply(1:nrow(dmrs), function(i) {
    row <- dmrs[i, ]
    # 該当するwindowを取得
    overlapping_windows <- win[chr == row$chr &
                                direction == row$direction &
                                start >= row$start &
                                end <= row$end]
    overlapping_windows <- overlapping_windows[order(start)]
    
    if (nrow(overlapping_windows) == 0) return(NULL)
    
    # 前方から弱いウィンドウを除去
    start_idx <- 1
    while (start_idx <= nrow(overlapping_windows) && 
           overlapping_windows$p[start_idx] > p_extend) {
      start_idx <- start_idx + 1
    }
    
    # 後方から弱いウィンドウを除去
    end_idx <- nrow(overlapping_windows)
    while (end_idx >= start_idx && 
           overlapping_windows$p[end_idx] > p_extend) {
      end_idx <- end_idx - 1
    }
    
    # 有効なウィンドウが残っているか
    if (start_idx > end_idx) return(NULL)
    
    trimmed_windows <- overlapping_windows[start_idx:end_idx]
    
    # 新しい境界を設定
    row$start <- min(trimmed_windows$start)
    row$end <- max(trimmed_windows$end)
    row$n_windows <- nrow(trimmed_windows)
    
    return(row)
  })
  
  trimmed <- data.table::rbindlist(trimmed[!sapply(trimmed, is.null)])
  message("[INFO] Trimmed: ", nrow(dmrs), " -> ", nrow(trimmed), " DMRs")
  return(trimmed)
}

# ポストフィルタリング: FP削減のための品質チェック
post_filter_dmrs <- function(dmrs, win, min_median_p, min_consistent_frac, p_seed) {
  if (nrow(dmrs) == 0) return(dmrs)
  
  message("[INFO] Applying post-filter: median_p <= ", min_median_p, 
          ", consistent_frac >= ", min_consistent_frac)
  
  keep <- sapply(1:nrow(dmrs), function(i) {
    row <- dmrs[i, ]
    # 該当するwindowを取得
    overlapping_windows <- win[chr == row$chr &
                                direction == row$direction &
                                start >= row$start &
                                end <= row$end]
    
    if (nrow(overlapping_windows) == 0) return(FALSE)
    
    # 品質チェック1: p値の中央値
    median_p <- median(overlapping_windows$p, na.rm = TRUE)
    
    # 品質チェック2: p < p_seedの窓の割合
    consistent_frac <- mean(overlapping_windows$p <= p_seed, na.rm = TRUE)
    
    # 両方の条件を満たすか
    median_p <= min_median_p && consistent_frac >= min_consistent_frac
  })
  
  filtered <- dmrs[keep, ]
  message("[INFO] Post-filter: ", nrow(dmrs), " -> ", nrow(filtered), 
          " DMRs (removed ", nrow(dmrs) - nrow(filtered), ")")
  
  filtered
}

merge_overlapping_dmrs <- function(dmrs, win, max_gap_bp = 0, p_seed = 0.05) {
  if (nrow(dmrs) == 0) return(dmrs)

  dmrs <- dmrs[order(chr, start, end)]
  merged <- list()
  i <- 1
  while (i <= nrow(dmrs)) {
    current <- dmrs[i, ]
    j <- i + 1
    while (j <= nrow(dmrs) &&
           dmrs$chr[j] == current$chr &&
           dmrs$direction[j] == current$direction) {
      gap <- dmrs$start[j] - current$end
      if (gap <= max_gap_bp) {
        current$end <- max(current$end, dmrs$end[j])
        j <- j + 1
      } else {
        break
      }
    }

    merged_windows <- win[chr == current$chr &
                           direction == current$direction &
                           start >= current$start &
                           end <= current$end]
    n_windows <- nrow(merged_windows)
    if (n_windows == 0) {
      i <- j
      next
    }

    combined_p <- p_stouffer(merged_windows$p)
    out_row <- data.table(
      chr = current$chr,
      start = current$start,
      end = current$end,
      n_windows = n_windows,
      direction = current$direction,
      combined_p = combined_p
    )

    if ("strong_frac" %in% names(dmrs)) {
      out_row[, strong_frac := mean(merged_windows$p <= p_seed, na.rm = TRUE)]
    }

    merged[[length(merged) + 1]] <- out_row
    i <- j
  }

  data.table::rbindlist(merged)
}

# Simes法でDMRを検出(双方向拡張)
detect_dmrs_simes <- function(win, p_seed, max_gap_bp, min_windows) {
  res <- list()
  i <- 1L
  n <- nrow(win)
  while (i <= n) {
    if (win$p[i] > p_seed) {
      i <- i + 1
      next
    }
    chr <- win$chr[i]
    direction <- win$direction[i]
    run_idx <- i
    
    # 前方拡張
    j <- i + 1
    while (j <= n && win$chr[j] == chr && win$direction[j] == direction) {
      gap <- win$start[j] - win$end[j - 1]
      if (gap > max_gap_bp) break
      run_idx <- c(run_idx, j)
      j <- j + 1
    }
    
    # 後方拡張
    k <- i - 1
    while (k >= 1 && win$chr[k] == chr && win$direction[k] == direction) {
      gap <- win$start[i] - win$end[k]
      if (gap > max_gap_bp) break
      run_idx <- c(k, run_idx)
      k <- k - 1
    }
    
    # Simes統合p値を計算
    simes_p <- p_simes(win$p[run_idx])
    if (length(run_idx) >= min_windows && simes_p <= p_seed) {
      res[[length(res) + 1]] <- list(
        chr = chr,
        start = min(win$start[run_idx]),
        end = max(win$end[run_idx]),
        n_windows = length(run_idx),
        direction = direction,
        combined_p = simes_p
      )
    }
    i <- max(run_idx) + 1
  }
  return(data.table::rbindlist(res))
}

# Stouffer法でDMRを検出(双方向拡張)
detect_dmrs_stouffer <- function(win, p_seed, p_extend, max_gap_bp, min_windows) {
  res <- list()
  i <- 1L
  n <- nrow(win)
  while (i <= n) {
    if (win$p[i] > p_seed) {
      i <- i + 1
      next
    }
    chr <- win$chr[i]
    direction <- win$direction[i]
    run_idx <- i
    
    # 前方拡張
    j <- i + 1
    while (j <= n && win$chr[j] == chr && win$direction[j] == direction) {
      gap <- win$start[j] - win$end[max(run_idx)]
      if (gap > max_gap_bp) break
      test_idx <- c(run_idx, j)
      z_scores <- qnorm(1 - win$p[test_idx])
      combined_z <- sum(z_scores) / sqrt(length(z_scores))
      stouffer_p <- 2 * pnorm(-abs(combined_z))
      if (stouffer_p <= p_extend) {
        run_idx <- test_idx
        j <- j + 1
      } else {
        break
      }
    }
    
    # 後方拡張
    k <- i - 1
    while (k >= 1 && win$chr[k] == chr && win$direction[k] == direction) {
      gap <- win$start[min(run_idx)] - win$end[k]
      if (gap > max_gap_bp) break
      test_idx <- c(k, run_idx)
      z_scores <- qnorm(1 - win$p[test_idx])
      combined_z <- sum(z_scores) / sqrt(length(z_scores))
      stouffer_p <- 2 * pnorm(-abs(combined_z))
      if (stouffer_p <= p_extend) {
        run_idx <- test_idx
        k <- k - 1
      } else {
        break
      }
    }
    
    # 最終的な統合p値を計算
    z_scores <- qnorm(1 - win$p[run_idx])
    combined_z <- sum(z_scores) / sqrt(length(z_scores))
    stouffer_p <- 2 * pnorm(-abs(combined_z))
    if (length(run_idx) >= min_windows && stouffer_p <= p_extend) {
      res[[length(res) + 1]] <- list(
        chr = chr,
        start = min(win$start[run_idx]),
        end = max(win$end[run_idx]),
        n_windows = length(run_idx),
        direction = direction,
        combined_p = stouffer_p
      )
    }
    i <- max(run_idx) + 1
  }
  return(data.table::rbindlist(res))
}

# Single Seed法でDMRを検出(単一窓シードからStouffer法で双方向拡張)
detect_dmrs_single_seed <- function(win, p_seed, p_extend, max_gap_bp, min_windows, min_delta = 0, 
                                    max_p_degradation = 1.2, max_final_p = 1.0, min_strong_windows = 0.5) {
  res <- list()
  i <- 1L
  n <- nrow(win)
  
  while (i <= n) {
    # 単一窓がp_seed以下ならシードとする
    if (win$p[i] > p_seed) {
      i <- i + 1
      next
    }
    
    chr <- win$chr[i]
    direction <- win$direction[i]
    run_idx <- i
    
    # 初期のStouffer p値（シード単体）
    current_stouffer_p <- win$p[i]
    
    # Stouffer法で前方拡張（p値悪化チェック + 効果量フィルタ）
    j <- i + 1
    while (j <= n && win$chr[j] == chr && win$direction[j] == direction) {
      gap <- win$start[j] - win$end[max(run_idx)]
      if (gap > max_gap_bp) break
      
      # 効果量フィルタ（設定されている場合は小さいdeltaをスキップ）
      if (min_delta > 0 && abs(win$delta[j]) < min_delta) {
        j <- j + 1
        next
      }
      
      test_idx <- c(run_idx, j)
      z_scores <- qnorm(1 - win$p[test_idx])
      combined_z <- sum(z_scores) / sqrt(length(z_scores))
      stouffer_p <- 2 * pnorm(-abs(combined_z))
      
      # p値悪化チェック: より厳格な閾値を使用
      if (stouffer_p <= p_extend && stouffer_p <= current_stouffer_p * max_p_degradation) {
        run_idx <- test_idx
        current_stouffer_p <- stouffer_p  # p値を更新
        j <- j + 1
      } else {
        break
      }
    }
    
    # Stouffer法で後方拡張（p値悪化チェック + 効果量フィルタ）
    k <- i - 1
    while (k >= 1 && win$chr[k] == chr && win$direction[k] == direction) {
      gap <- win$start[min(run_idx)] - win$end[k]
      if (gap > max_gap_bp) break
      
      # 効果量フィルタ
      if (min_delta > 0 && abs(win$delta[k]) < min_delta) {
        k <- k - 1
        next
      }
      
      test_idx <- c(k, run_idx)
      z_scores <- qnorm(1 - win$p[test_idx])
      combined_z <- sum(z_scores) / sqrt(length(z_scores))
      stouffer_p <- 2 * pnorm(-abs(combined_z))
      
      # p値悪化チェック: より厳格な閾値を使用
      if (stouffer_p <= p_extend && stouffer_p <= current_stouffer_p * max_p_degradation) {
        run_idx <- test_idx
        current_stouffer_p <- stouffer_p
        k <- k - 1
      } else {
        break
      }
    }
    
    # 最終的な統合p値を計算
    z_scores_final <- qnorm(1 - win$p[run_idx])
    combined_z_final <- sum(z_scores_final) / sqrt(length(z_scores_final))
    final_stouffer_p <- 2 * pnorm(-abs(combined_z_final))
    
    # 品質チェック: 強いシグナルのウィンドウ割合
    strong_windows_frac <- sum(win$p[run_idx] <= p_seed) / length(run_idx)
    
    # DMRを追加（フィルタリング条件を適用）
    if (length(run_idx) >= min_windows && 
        final_stouffer_p <= max_final_p && 
        strong_windows_frac >= min_strong_windows) {
      res[[length(res) + 1]] <- list(
        chr = chr,
        start = min(win$start[run_idx]),
        end = max(win$end[run_idx]),
        n_windows = length(run_idx),
        direction = direction,
        combined_p = final_stouffer_p,
        strong_frac = strong_windows_frac
      )
    }
    
    i <- max(run_idx) + 1
  }
  
  return(data.table::rbindlist(res))
}

# Stouffer_multi_seed法: 複数ウィンドウでStoufferシードを作成（truth DMR検出率向上）
detect_dmrs_stouffer_multi_seed <- function(win, p_seed, p_extend, max_gap_bp, min_windows, min_delta = 0, 
                                             max_p_degradation = 1.2, max_final_p = 1.0, min_strong_windows = 0.5, 
                                             seed_min_windows = 2) {
  res <- list()
  i <- 1L
  n <- nrow(win)
  
  while (i <= n) {
    chr <- win$chr[i]
    direction <- win$direction[i]
    
    # 複数ウィンドウでStoufferシードを作成
    seed_idx <- i
    j <- i + 1
    while (j <= n && win$chr[j] == chr && win$direction[j] == direction) {
      gap <- win$start[j] - win$end[j - 1]
      if (gap > max_gap_bp) break
      seed_idx <- c(seed_idx, j)
      j <- j + 1
      
      # シードが十分なサイズになったらStouffer検定
      if (length(seed_idx) >= seed_min_windows) {
        z_scores <- qnorm(1 - win$p[seed_idx])
        combined_z <- sum(z_scores) / sqrt(length(z_scores))
        stouffer_seed_p <- 2 * pnorm(-abs(combined_z))
        
        if (stouffer_seed_p <= p_seed) {
          # 有効なシードとして採用
          break
        }
      }
    }
    
    # シードが条件を満たさない場合は次へ
    if (length(seed_idx) < seed_min_windows) {
      i <- max(seed_idx) + 1
      next
    }
    
    z_scores_seed <- qnorm(1 - win$p[seed_idx])
    combined_z_seed <- sum(z_scores_seed) / sqrt(length(z_scores_seed))
    stouffer_seed_p <- 2 * pnorm(-abs(combined_z_seed))
    
    if (stouffer_seed_p > p_seed) {
      i <- max(seed_idx) + 1
      next
    }
    
    # シードから前方・後方拡張
    run_idx <- seed_idx
    current_stouffer_p <- stouffer_seed_p
    
    # 前方拡張
    j <- max(run_idx) + 1
    while (j <= n && win$chr[j] == chr && win$direction[j] == direction) {
      gap <- win$start[j] - win$end[max(run_idx)]
      if (gap > max_gap_bp) break
      
      if (min_delta > 0 && abs(win$delta[j]) < min_delta) {
        j <- j + 1
        next
      }
      
      test_idx <- c(run_idx, j)
      z_scores <- qnorm(1 - win$p[test_idx])
      combined_z <- sum(z_scores) / sqrt(length(z_scores))
      stouffer_p <- 2 * pnorm(-abs(combined_z))
      
      if (stouffer_p <= p_extend && stouffer_p <= current_stouffer_p * max_p_degradation) {
        run_idx <- test_idx
        current_stouffer_p <- stouffer_p
        j <- j + 1
      } else {
        break
      }
    }
    
    # 後方拡張
    k <- min(seed_idx) - 1
    while (k >= 1 && win$chr[k] == chr && win$direction[k] == direction) {
      gap <- win$start[min(run_idx)] - win$end[k]
      if (gap > max_gap_bp) break
      
      if (min_delta > 0 && abs(win$delta[k]) < min_delta) {
        k <- k - 1
        next
      }
      
      test_idx <- c(k, run_idx)
      z_scores <- qnorm(1 - win$p[test_idx])
      combined_z <- sum(z_scores) / sqrt(length(z_scores))
      stouffer_p <- 2 * pnorm(-abs(combined_z))
      
      if (stouffer_p <= p_extend && stouffer_p <= current_stouffer_p * max_p_degradation) {
        run_idx <- test_idx
        current_stouffer_p <- stouffer_p
        k <- k - 1
      } else {
        break
      }
    }
    
    # 最終統合p値
    z_scores_final <- qnorm(1 - win$p[run_idx])
    combined_z_final <- sum(z_scores_final) / sqrt(length(z_scores_final))
    final_stouffer_p <- 2 * pnorm(-abs(combined_z_final))
    
    # 品質チェック
    strong_windows_frac <- sum(win$p[run_idx] <= p_seed) / length(run_idx)
    
    if (length(run_idx) >= min_windows && 
        final_stouffer_p <= max_final_p && 
        strong_windows_frac >= min_strong_windows) {
      res[[length(res) + 1]] <- list(
        chr = chr,
        start = min(win$start[run_idx]),
        end = max(win$end[run_idx]),
        n_windows = length(run_idx),
        direction = direction,
        combined_p = final_stouffer_p,
        strong_frac = strong_windows_frac
      )
    }
    
    i <- max(run_idx) + 1
  }
  
  return(data.table::rbindlist(res))
}

detect_dmrs_combined <- function(win, p_seed, p_extend, max_gap_bp, min_windows, min_delta = 0) {
  res <- list()
  i <- 1L
  n <- nrow(win)
  
  while (i <= n) {
    if (win$p[i] > p_seed) {
      i <- i + 1
      next
    }
    chr <- win$chr[i]
    direction <- win$direction[i]
    
    # Simes法でシード（起点）を検出
    seed_idx <- i
    j <- i + 1
    while (j <= n && win$chr[j] == chr && win$direction[j] == direction) {
      gap <- win$start[j] - win$end[j - 1]
      if (gap > max_gap_bp) break
      seed_idx <- c(seed_idx, j)
      j <- j + 1
    }
    
    # Simes統合p値で起点として適切か判定
    simes_seed_p <- p_simes(win$p[seed_idx])
    if (simes_seed_p > p_seed || length(seed_idx) < min_windows) {
      i <- max(seed_idx) + 1
      next
    }
    
    # 初期のStouffer p値を計算
    z_scores_init <- qnorm(1 - win$p[seed_idx])
    combined_z_init <- sum(z_scores_init) / sqrt(length(z_scores_init))
    current_stouffer_p <- 2 * pnorm(-abs(combined_z_init))
    
    # Stouffer法で前方拡張（p値悪化チェック + 効果量フィルタ）
    run_idx <- seed_idx
    j <- max(run_idx) + 1
    while (j <= n && win$chr[j] == chr && win$direction[j] == direction) {
      gap <- win$start[j] - win$end[max(run_idx)]
      if (gap > max_gap_bp) break
      
      # 効果量フィルタ（設定されている場合は小さいdeltaをスキップ）
      if (min_delta > 0 && abs(win$delta[j]) < min_delta) {
        j <- j + 1
        next
      }
      
      test_idx <- c(run_idx, j)
      z_scores <- qnorm(1 - win$p[test_idx])
      combined_z <- sum(z_scores) / sqrt(length(z_scores))
      stouffer_p <- 2 * pnorm(-abs(combined_z))
      
      # p値悪化チェック: 新しいp値が現在のp値より悪化（大きくなる）する場合は拡張停止
      if (stouffer_p <= p_extend && stouffer_p <= current_stouffer_p * 1.5) {
        run_idx <- test_idx
        current_stouffer_p <- stouffer_p  # p値を更新
        j <- j + 1
      } else {
        break
      }
    }
    
    # Stouffer法で後方拡張（p値悪化チェック + 効果量フィルタ）
    k <- min(seed_idx) - 1
    while (k >= 1 && win$chr[k] == chr && win$direction[k] == direction) {
      gap <- win$start[min(run_idx)] - win$end[k]
      if (gap > max_gap_bp) break
      
      # 効果量フィルタ
      if (min_delta > 0 && abs(win$delta[k]) < min_delta) {
        k <- k - 1
        next
      }
      
      test_idx <- c(k, run_idx)
      z_scores <- qnorm(1 - win$p[test_idx])
      combined_z <- sum(z_scores) / sqrt(length(z_scores))
      stouffer_p <- 2 * pnorm(-abs(combined_z))
      
      # p値悪化チェック
      if (stouffer_p <= p_extend && stouffer_p <= current_stouffer_p * 1.5) {
        run_idx <- test_idx
        current_stouffer_p <- stouffer_p
        k <- k - 1
      } else {
        break
      }
    }
    
    # 最終的な統合p値を計算
    final_simes_p <- p_simes(win$p[run_idx])
    z_scores_final <- qnorm(1 - win$p[run_idx])
    combined_z_final <- sum(z_scores_final) / sqrt(length(z_scores_final))
    final_stouffer_p <- 2 * pnorm(-abs(combined_z_final))
    
    # DMRを追加
    res[[length(res) + 1]] <- list(
      chr = chr,
      start = min(win$start[run_idx]),
      end = max(win$end[run_idx]),
      n_windows = length(run_idx),
      direction = direction,
      simes_p = final_simes_p,
      stouffer_p = final_stouffer_p
    )
    
    i <- max(run_idx) + 1
  }
  
  return(data.table::rbindlist(res))
}

#========================= Main Script ==========================================
win <- load_data(opt$windows)

valid_modes <- c("Simes", "Stouffer", "single_seed", "multi_seed", "hybrid_seed")
if (!opt$`merge-mode` %in% valid_modes) {
  stop(sprintf("Invalid --merge-mode '%s'. Must be one of: %s",
               opt$`merge-mode`, paste(valid_modes, collapse = ", ")))
}

if (opt$`merge-mode` == "Simes") {
  dmrs_simes <- detect_dmrs_simes(win, opt$`p-seed`, opt$`max-gap-bp`, opt$`min-windows`)
  dmrs_simes <- add_delta_summary(dmrs_simes, win)  # Delta値を追加
  fwrite(dmrs_simes, sprintf("%s_dmrs_simes.tsv", opt$`out-prefix`), sep = "\t")
  write_bed(dmrs_simes, sprintf("%s_dmrs_simes.bed", opt$`out-prefix`))  # BED出力
} else if (opt$`merge-mode` == "Stouffer") {
  dmrs_stouffer <- detect_dmrs_stouffer(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`)
  if (opt$`post-filter`) {
    dmrs_stouffer <- post_filter_dmrs(dmrs_stouffer, win, opt$`min-median-p`, 
                                      opt$`min-consistent-frac`, opt$`p-seed`)
  }
  dmrs_stouffer <- add_delta_summary(dmrs_stouffer, win)  # Delta値を追加
  fwrite(dmrs_stouffer, sprintf("%s_dmrs_stouffer.tsv", opt$`out-prefix`), sep = "\t")
  write_bed(dmrs_stouffer, sprintf("%s_dmrs_stouffer.bed", opt$`out-prefix`))  # BED出力
} else if (opt$`merge-mode` == "multi_seed") {
  # 適応的delta閾値を使用する場合: まず緩い条件でDMRを検出
  if (opt$`adaptive-delta`) {
    message("[INFO] Using adaptive delta threshold...")
    # 第1段階: min-delta=0でDMR候補を検出
    dmrs_initial <- detect_dmrs_stouffer_multi_seed(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`, 
                                                     min_delta = 0, opt$`max-p-degradation`, opt$`max-final-p`, opt$`min-strong-windows`,
                                                     seed_min_windows = opt$`seed-min-windows`)
    
    # 第2段階: DMR候補のdelta分布から適応的閾値を計算
    adaptive_min_delta <- calculate_adaptive_delta_threshold(dmrs_initial, win, 
                                                              method = opt$`adaptive-delta-method`, 
                                                              ratio = opt$`adaptive-delta-ratio`)
    
    # 第3段階: 適応的閾値で再検出
    dmrs_stouffer_ms <- detect_dmrs_stouffer_multi_seed(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`, 
                                                         adaptive_min_delta, opt$`max-p-degradation`, opt$`max-final-p`, opt$`min-strong-windows`,
                                                         seed_min_windows = opt$`seed-min-windows`)
  } else {
    dmrs_stouffer_ms <- detect_dmrs_stouffer_multi_seed(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`, 
                                                         opt$`min-delta`, opt$`max-p-degradation`, opt$`max-final-p`, opt$`min-strong-windows`,
                                                         seed_min_windows = opt$`seed-min-windows`)
  }
  
  if (opt$`trim-weak-edges`) {
    dmrs_stouffer_ms <- trim_dmr_edges(dmrs_stouffer_ms, win, opt$`p-extend`)
  }
  
  if (opt$`min-dmr-length` > 0) {
    before_n <- nrow(dmrs_stouffer_ms)
    dmrs_stouffer_ms <- dmrs_stouffer_ms[end - start >= opt$`min-dmr-length`, ]
    message("[INFO] DMR length filter: ", before_n, " -> ", nrow(dmrs_stouffer_ms))
  }
  
  if (opt$`max-median-p` < 1.0 && nrow(dmrs_stouffer_ms) > 0) {
    before_n <- nrow(dmrs_stouffer_ms)
    dmrs_stouffer_ms <- lapply(1:nrow(dmrs_stouffer_ms), function(i) {
      row <- dmrs_stouffer_ms[i, ]
      overlapping_windows <- win[chr == row$chr & direction == row$direction &
                                   start >= row$start & end <= row$end]
      median_p <- median(overlapping_windows$p, na.rm = TRUE)
      if (median_p <= opt$`max-median-p`) return(row) else return(NULL)
    })
    dmrs_stouffer_ms <- data.table::rbindlist(dmrs_stouffer_ms[!sapply(dmrs_stouffer_ms, is.null)])
    message("[INFO] Median p-value filter: ", before_n, " -> ", nrow(dmrs_stouffer_ms))
  }
  
  if (opt$`post-filter`) {
    dmrs_stouffer_ms <- post_filter_dmrs(dmrs_stouffer_ms, win, opt$`min-median-p`, 
                                         opt$`min-consistent-frac`, opt$`p-seed`)
  }
  if (opt$`merge-overlaps`) {
    dmrs_stouffer_ms <- merge_overlapping_dmrs(
      dmrs_stouffer_ms,
      win,
      max_gap_bp = opt$`merge-overlaps-gap`,
      p_seed = opt$`p-seed`
    )
  }
  dmrs_multi_seed <- add_delta_summary(dmrs_stouffer_ms, win)
  fwrite(dmrs_multi_seed, sprintf("%s_dmrs_multi_seed.tsv", opt$`out-prefix`), sep = "\t")
  write_bed(dmrs_multi_seed, sprintf("%s_dmrs_multi_seed.bed", opt$`out-prefix`))
} else if (opt$`merge-mode` == "hybrid_seed") {
  dmrs_hybrid_seed <- detect_dmrs_combined(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`, opt$`min-delta`)
  if (opt$`post-filter`) {
    dmrs_hybrid_seed <- post_filter_dmrs(dmrs_hybrid_seed, win, opt$`min-median-p`,
                                         opt$`min-consistent-frac`, opt$`p-seed`)
  }
  dmrs_hybrid_seed <- add_delta_summary(dmrs_hybrid_seed, win)  # Delta値を追加
  fwrite(dmrs_hybrid_seed, sprintf("%s_dmrs_hybrid_seed.tsv", opt$`out-prefix`), sep = "\t")
  write_bed(dmrs_hybrid_seed, sprintf("%s_dmrs_hybrid_seed.bed", opt$`out-prefix`))  # BED出力
} else if (opt$`merge-mode` == "single_seed") {
  dmrs_single_seed <- detect_dmrs_single_seed(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`, 
                                              opt$`min-delta`, opt$`max-p-degradation`, opt$`max-final-p`, opt$`min-strong-windows`)
  
  # 境界トリミング（オプション）
  if (opt$`trim-weak-edges`) {
    dmrs_single_seed <- trim_dmr_edges(dmrs_single_seed, win, opt$`p-extend`)
  }
  
  # DMRサイズフィルタ
  if (opt$`min-dmr-length` > 0) {
    before_n <- nrow(dmrs_single_seed)
    dmrs_single_seed <- dmrs_single_seed[end - start >= opt$`min-dmr-length`, ]
    message("[INFO] DMR length filter (>= ", opt$`min-dmr-length`, " bp): ", before_n, " -> ", nrow(dmrs_single_seed))
  }
  
  # 中央値p値フィルタ
  if (opt$`max-median-p` < 1.0) {
    before_n <- nrow(dmrs_single_seed)
    dmrs_single_seed <- lapply(1:nrow(dmrs_single_seed), function(i) {
      row <- dmrs_single_seed[i, ]
      overlapping_windows <- win[chr == row$chr & direction == row$direction &
                                   start >= row$start & end <= row$end]
      median_p <- median(overlapping_windows$p, na.rm = TRUE)
      if (median_p <= opt$`max-median-p`) return(row) else return(NULL)
    })
    dmrs_single_seed <- data.table::rbindlist(dmrs_single_seed[!sapply(dmrs_single_seed, is.null)])
    message("[INFO] Median p-value filter (<= ", opt$`max-median-p`, "): ", before_n, " -> ", nrow(dmrs_single_seed))
  }
  
  if (opt$`post-filter`) {
    dmrs_single_seed <- post_filter_dmrs(dmrs_single_seed, win, opt$`min-median-p`, 
                                         opt$`min-consistent-frac`, opt$`p-seed`)
  }
  dmrs_single_seed <- add_delta_summary(dmrs_single_seed, win)  # Delta値を追加
  fwrite(dmrs_single_seed, sprintf("%s_dmrs_single_seed.tsv", opt$`out-prefix`), sep = "\t")
  write_bed(dmrs_single_seed, sprintf("%s_dmrs_single_seed.bed", opt$`out-prefix`))  # BED出力
}

message("[INFO] DMR analysis complete.")
