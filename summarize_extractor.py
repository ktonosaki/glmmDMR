#!/usr/bin/env python3

import argparse
import gzip
import os
import io
from pathlib import Path
import pandas as pd
from tqdm import tqdm
import multiprocessing
from functools import partial

def read_methylation_file(file, context):
    with gzip.open(file, 'rt') as f:
        filtered_lines = [line for line in f if not line.startswith("Bismark methylation extractor version")]

    df = pd.read_csv(
        io.StringIO(''.join(filtered_lines)),
        sep='\t',
        header=None,
        names=["read_id", "strand", "chr", "pos", "code"],
        usecols=["chr", "pos", "strand", "code"],
        dtype={"chr": str, "pos": int, "strand": str, "code": str}
    )

#    df['meth'] = df['code'].isin(['Z', 'X', 'H']).astype(int)
#    df['unmeth'] = df['code'].isin(['z', 'x', 'h']).astype(int)

    # codeごとに1レコード1リードとみなし、各codeの出現数をカウントする
    df['meth'] = df['code'].isin(['Z', 'X', 'H'])  # メチル化コード
    df['unmeth'] = df['code'].isin(['z', 'x', 'h'])  # 非メチル化コード

    # bool → int に変換
    df['meth'] = df['meth'].astype(int)
    df['unmeth'] = df['unmeth'].astype(int)

    df = df.groupby(['chr', 'pos', 'strand'], as_index=False)[['meth', 'unmeth']].sum()
    df['context'] = context
    df = df[['chr', 'pos', 'strand', 'meth', 'unmeth', 'context']]
    return df

def process_context(context, files, output_prefix, threads, keep_temp):
#    dfs = []
    print(f"Processing context {context} with {len(files)} files...")
#    with multiprocessing.Pool(threads) as pool:
#        dfs = list(tqdm(pool.imap(partial(read_methylation_file, context=context), files), total=len(files)))
    with multiprocessing.Pool(threads) as pool:
        dfs = list(tqdm(pool.imap(partial(read_methylation_file, context=context), files),
                        total=len(files), desc=f"{context}", dynamic_ncols=True))

    merged_df = pd.concat(dfs, axis=0)
    merged_df = merged_df.groupby(['chr', 'pos', 'strand', 'context'], as_index=False)[['meth', 'unmeth']].sum()
    tmp_file = f"{output_prefix}.{context}.tmp.tsv"
    merged_df.to_csv(tmp_file, sep='\t', index=False)

#    if not keep_temp:
#        for f in files:
#            os.remove(f)
    # 注意：ここで files を削除するのは危険。通常 keep_temp による削除は tmp_file に適用する方が適切
    if not keep_temp:
        os.remove(tmp_file)

    return tmp_file

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input directory containing bismark_methylation_extractor OUTPUT files")
    parser.add_argument("-o", "--output", required=True, help="Output summarized tsv.gz file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use (default: 4)")
    parser.add_argument("--keep-temp", action="store_true", help="Keep intermediate temp files")
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_file = args.output
    threads = args.threads
    keep_temp = args.keep_temp

    # 固定の順序
    context_order = ['CpG', 'CHG', 'CHH']
    context_files = {ctx: [] for ctx in context_order}

    print("Scanning files...")
    for file in sorted(input_dir.glob("*.txt.gz")):
        for ctx in context_order:
            if file.name.startswith(ctx + "_"):
                context_files[ctx].append(str(file))

    temp_files = []
    for ctx in context_order:
        if context_files[ctx]:
            temp_file = process_context(ctx, context_files[ctx], output_file, threads, keep_temp)
            temp_files.append(temp_file)

    print("Combining all contexts...")
    combined_df = pd.concat([pd.read_csv(f, sep='\t') for f in temp_files], axis=0)
    combined_df = combined_df[['chr', 'pos', 'strand', 'meth', 'unmeth', 'context']]
    combined_df.to_csv(output_file, sep='\t', index=False, compression='gzip')

    print("\n=== Strand-wise summary ===")
    print(combined_df.groupby('strand')[['meth', 'unmeth']].sum())

if __name__ == "__main__":
    main()
