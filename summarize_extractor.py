#!/usr/bin/env python3

import argparse
import gzip
import os
import multiprocessing
from functools import partial
from pathlib import Path
import pandas as pd
import io
from tqdm import tqdm

def read_methylation_file(file, context):
    file = Path(file)
    suffix = file.name.split("_")[1]  # e.g. OT, OB, CTOT, CTOB

    if suffix in ['OT', 'CTOT']:
        true_strand = '+'
    elif suffix in ['OB', 'CTOB']:
        true_strand = '-'
    else:
        raise ValueError(f"Unrecognized strand suffix in filename: {file.name}")

    with gzip.open(file, 'rt') as f:
        # Skip all Bismark header lines (more robust)
        filtered_lines = [line for line in f if not line.startswith("Bismark")]
        
    if not filtered_lines:
        print(f"Warning: Empty file {file.name}")
        return pd.DataFrame(columns=["chr", "pos", "strand", "meth", "unmeth", "context"])
    
    df = pd.read_csv(
        io.StringIO(''.join(filtered_lines)),
        sep='\t',
        header=None,
        names=["read_id", "strand", "chr", "pos", "code"],
        usecols=["chr", "pos", "code"],
        dtype={"chr": str, "pos": int, "code": str}
    )

    # Set strand based on filename
    df['strand'] = true_strand

    df['meth'] = df['code'].isin(['Z', 'X', 'H'])
    df['unmeth'] = df['code'].isin(['z', 'x', 'h'])

    df = df.groupby(['chr', 'pos', 'strand'], as_index=False)[['meth', 'unmeth']].sum()
    df['context'] = context
    df = df[['chr', 'pos', 'strand', 'meth', 'unmeth', 'context']]
    return df


def process_context(context, files, output_prefix, threads, keep_temp):
    print(f"Processing context {context} with {len(files)} files...")
    with multiprocessing.Pool(threads) as pool:
        dfs = list(tqdm(pool.imap(partial(read_methylation_file, context=context), files), total=len(files)))
    
    # Filter out empty DataFrames
    dfs = [df for df in dfs if not df.empty]
    if not dfs:
        print(f"Warning: No data for context {context}")
        return None
    
    merged_df = pd.concat(dfs, axis=0)
    merged_df = merged_df.groupby(['chr', 'pos', 'strand', 'context'], as_index=False)[['meth', 'unmeth']].sum()
    
    tmp_file = f"{output_prefix}.{context}.tmp.tsv"
    merged_df.to_csv(tmp_file, sep='\t', index=False)
    
    return tmp_file


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input directory containing deduplicated.txt.gz files")
    parser.add_argument("-o", "--output", required=True, help="Output summarized tsv.gz file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use (default: 4)")
    parser.add_argument("--keep-temp", action="store_true", help="Keep intermediate temp files")
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_file = args.output
    threads = args.threads
    keep_temp = args.keep_temp

    # Validate input directory
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    
    # Create output directory if needed
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    context_order = ['CpG', 'CHG', 'CHH']
    context_files = {ctx: [] for ctx in context_order}

    print("Scanning files...")
    for file in sorted(input_dir.glob("*.txt.gz")):
        for ctx in context_order:
            if file.name.startswith(ctx + "_"):
                context_files[ctx].append(str(file))
    
    # Check if any files found
    total_files = sum(len(files) for files in context_files.values())
    if total_files == 0:
        raise FileNotFoundError(f"No .txt.gz files found in {input_dir}")

    temp_files = []
    for ctx in context_order:
        if context_files[ctx]:
            temp_file = process_context(ctx, context_files[ctx], output_file, threads, keep_temp)
            if temp_file:
                temp_files.append(temp_file)
    
    if not temp_files:
        raise ValueError("No data processed from any context")

    print("Combining all contexts...")
    combined_df = pd.concat(
        [pd.read_csv(f, sep='\t', dtype={'chr': str}) for f in temp_files],
        axis=0
    )
    combined_df = combined_df[['chr', 'pos', 'strand', 'meth', 'unmeth', 'context']]
    combined_df.to_csv(output_file, sep='\t', index=False, compression='gzip')
    
    # Cleanup temp files unless --keep-temp specified
    if not keep_temp:
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)
                print(f"Removed temp file: {temp_file}")

    print(f"\nOutput written to: {output_file}")
    print("\n=== Strand-wise summary ===")
    print(combined_df.groupby('strand')[['meth', 'unmeth']].sum())

if __name__ == "__main__":
    main()
