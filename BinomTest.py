#!/usr/bin/env python3

import argparse
import gzip
import pandas as pd
import numpy as np
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests
import multiprocessing as mp
from functools import partial
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

def binom_test_row(row, null_prob):
    n_total = row['meth'] + row['unmeth']
    return binomtest(row['meth'], n=n_total, p=null_prob, alternative='two-sided').pvalue

def process_chunk(df_chunk, null_prob):
    df_chunk['pval'] = df_chunk.apply(binom_test_row, axis=1, null_prob=null_prob)
    return df_chunk

def main():
    parser = argparse.ArgumentParser(description="Binomial test for methylation sites with FDR correction")
    parser.add_argument("-i", "--input", required=True, help="Input summarized file (tsv.gz)")
    parser.add_argument("-o", "--output", required=True, help="Output filtered file (tsv.gz)")
    parser.add_argument("--null_prob", type=float, default=None, help="Null hypothesis probability (default: estimated from nonconv_chr or 0.5)")
    parser.add_argument("--fdr_threshold", type=float, default=0.05, help="FDR threshold (default: 0.05)")
    parser.add_argument("--min_coverage", type=int, default=0, help="Minimum coverage to retain site (default: 0)")
    parser.add_argument("--nonconv_chr", type=str, default=None, help="Non-conversion control chromosome (optional, e.g., 'chloroplast')")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads (default: 4)")
    args = parser.parse_args()

    # Validate input
    from pathlib import Path
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")
    
    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Reading input: {args.input}")
    df = pd.read_csv(args.input, sep='\t', compression='infer', low_memory=False, dtype={'chr': str})
    
    if df.empty:
        raise ValueError("Input file is empty")
    
    print(f"Loaded {len(df):,} sites")
    
    # Ensure required columns
    if 'start' not in df.columns and 'pos' in df.columns:
        df['start'] = df['pos']
    df['end'] = df['start'] + 1
    df['coverage'] = df['meth'] + df['unmeth']
    # Filter by coverage (single pass)
    min_cov = max(1, args.min_coverage)  # Ensure at least 1 to avoid division by zero
    df = df[df['coverage'] >= min_cov]
    print(f"After coverage filter (>={min_cov}): {len(df):,} sites")
    
    if df.empty:
        raise ValueError(f"No sites remaining after coverage filter (min={min_cov})")

    # Estimate or use null probability
    null_prob = args.null_prob
    if null_prob is None and args.nonconv_chr:
        subset = df[df['chr'] == args.nonconv_chr]
        if not subset.empty:
            null_prob = subset['meth'].sum() / (subset['meth'].sum() + subset['unmeth'].sum())
            print(f"Estimated null_prob from {args.nonconv_chr}: {null_prob:.4f}")
        else:
            print(f"Warning: No data found for nonconv_chr '{args.nonconv_chr}', using default 0.5")
            null_prob = 0.5
    elif null_prob is None:
        null_prob = 0.5
        print(f"Using default null_prob: {null_prob}")
    else:
        print(f"Using specified null_prob: {null_prob}")

    # Parallel binomial test
    print(f"Running binomial tests (null_prob={null_prob:.4f}, threads={args.threads})...")
    chunks = np.array_split(df, args.threads)
    with mp.Pool(args.threads) as pool:
        results = list(tqdm(pool.imap(partial(process_chunk, null_prob=null_prob), chunks), 
                           total=len(chunks), desc="Processing chunks"))

    df = pd.concat(results, axis=0)
    
    # FDR correction
    print("Applying FDR correction...")
    df['FDR'] = multipletests(df['pval'], method='fdr_bh')[1]
    
    n_sig = (df['FDR'] <= args.fdr_threshold).sum()
    print(f"Significant sites (FDR <= {args.fdr_threshold}): {n_sig:,} / {len(df):,} ({100*n_sig/len(df):.2f}%)")
    
    # Set meth=0 for non-significant sites (Weighted methylation level approach)
    df.loc[df['FDR'] > args.fdr_threshold, 'meth'] = 0
    print(f"Non-significant sites set to meth=0 for Weighted methylation level calculation")
    
    # Output
    df = df[['chr', 'pos', 'strand', 'meth', 'unmeth', 'context']]
    df.to_csv(args.output, sep='\t', index=False, compression='gzip')
    print(f"\nOutput written to: {args.output}")

if __name__ == "__main__":
    main()
