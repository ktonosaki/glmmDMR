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
    return binomtest(int(row['meth']), n=int(row['meth'] + row['unmeth']), p=null_prob, alternative='two-sided').pvalue

def process_chunk(df_chunk, null_prob):
    df_chunk['pval'] = df_chunk.apply(binom_test_row, axis=1, null_prob=null_prob)
    return df_chunk

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input summarized file (tsv.gz)")
    parser.add_argument("-o", "--output", required=True, help="Output filtered file (tsv.gz)")
    parser.add_argument("--null_prob", type=float, default=None, help="Null hypothesis probability (default: estimated from nonconv_chr)")
    parser.add_argument("--fdr_threshold", type=float, default=0.05, help="FDR threshold (default: 0.05)")
    parser.add_argument("--min_coverage", type=int, default=0, help="Minimum coverage to retain site (default: 0)")
    parser.add_argument("--nonconv_chr", type=str, default=None, help="Non-conversion control chromosome (optional)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads (default: 4)")
    args = parser.parse_args()

#    df = pd.read_csv(args.input, sep='\t', compression='infer')
    df = pd.read_csv(args.input, sep='\t', compression='infer', low_memory=False)
    if 'start' not in df.columns and 'pos' in df.columns:
        df['start'] = df['pos']
    df['end'] = df['start'] + 1
    df['coverage'] = df['meth'] + df['unmeth']
    df = df[df['coverage'] > 0]

    if args.min_coverage > 0:
        df = df[df['coverage'] >= args.min_coverage]

    null_prob = args.null_prob
    if null_prob is None and args.nonconv_chr:
        subset = df[df['chr'] == args.nonconv_chr]
        if not subset.empty:
            null_prob = subset['meth'].sum() / (subset['meth'].sum() + subset['unmeth'].sum())
        else:
            raise ValueError("No data found for specified nonconv_chr")
    elif null_prob is None:
        null_prob = 0.5

    chunks = np.array_split(df, args.threads)
    with mp.Pool(args.threads) as pool:
        results = list(tqdm(pool.imap(partial(process_chunk, null_prob=null_prob), chunks), total=len(chunks)))

    df = pd.concat(results, axis=0)
    df['FDR'] = multipletests(df['pval'], method='fdr_bh')[1]
    df.loc[df['FDR'] > args.fdr_threshold, 'meth'] = 0
    df = df[['chr', 'pos', 'strand', 'meth', 'unmeth', 'context']]
    df.to_csv(args.output, sep='\t', index=False, compression='gzip')

if __name__ == "__main__":
    main()
