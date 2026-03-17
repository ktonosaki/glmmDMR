#!/usr/bin/env python3

import argparse
import math
import sys
from typing import List

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None

try:
    import pyBigWig
except ImportError as exc:  # pragma: no cover
    raise SystemExit("pyBigWig is required. Install with: pip install pyBigWig") from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute per-bin variance across multiple BigWig files and write a BigWig."
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Input BigWig files (space-separated).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output BigWig path.",
    )
    parser.add_argument(
        "--bin-size",
        type=int,
        default=200,
        help="Bin size in bp [default: 200].",
    )
    parser.add_argument(
        "--min-tracks",
        type=int,
        default=2,
        help="Minimum number of tracks required to compute variance [default: 2].",
    )
    parser.add_argument(
        "--norm",
        choices=["none", "log2p1"],
        default="none",
        help="Normalization before variance: none or log2p1 [default: none].",
    )
    return parser.parse_args()


def normalize(values, mode: str):
    if mode == "log2p1":
        if np is not None:
            return np.log2(values + 1.0)
        return [math.log2(v + 1.0) for v in values]
    return values


def stats_per_bin(bw, chrom: str, chrom_len: int, bin_size: int) -> List[float]:
    bins = int(math.ceil(chrom_len / bin_size))
    values = bw.stats(chrom, 0, chrom_len, type="mean", nBins=bins)
    return values


def main() -> int:
    args = parse_args()

    if args.bin_size <= 0:
        raise SystemExit("--bin-size must be > 0")

    if args.min_tracks <= 0:
        raise SystemExit("--min-tracks must be > 0")

    if len(args.inputs) < 2:
        sys.stderr.write("[WARN] Only one input provided; variance will be 0 or NA.\n")

    sys.stderr.write(f"[INFO] Inputs: {len(args.inputs)} files\n")
    sys.stderr.write(f"[INFO] Bin size: {args.bin_size} bp\n")
    sys.stderr.write(f"[INFO] Min tracks: {args.min_tracks}\n")
    sys.stderr.write(f"[INFO] Norm: {args.norm}\n")

    bws = [pyBigWig.open(path) for path in args.inputs]
    try:
        chroms = bws[0].chroms()
        bw_out = pyBigWig.open(args.output, "w")
        bw_out.addHeader(list(chroms.items()))

        for chrom, chrom_len in chroms.items():
            if chrom_len <= 0:
                continue

            per_track = []
            for bw in bws:
                values = stats_per_bin(bw, chrom, chrom_len, args.bin_size)
                per_track.append(values)

            bins = len(per_track[0]) if per_track else 0
            if bins == 0:
                continue

            if np is None:
                sys.stderr.write("[ERROR] numpy is required for variance calculation.\n")
                return 1

            arr = np.array(per_track, dtype=float)
            arr = normalize(arr, args.norm)

            counts = np.sum(~np.isnan(arr), axis=0)
            var = np.nanvar(arr, axis=0, ddof=0)
            var[counts < args.min_tracks] = np.nan

            starts = []
            ends = []
            values_out = []
            for i in range(bins):
                v = var[i]
                if np.isnan(v):
                    continue
                start = i * args.bin_size
                end = min((i + 1) * args.bin_size, chrom_len)
                starts.append(start)
                ends.append(end)
                values_out.append(float(v))

            if starts:
                bw_out.addEntries(
                    [chrom] * len(starts),
                    starts,
                    ends=ends,
                    values=values_out,
                )

        bw_out.close()
    finally:
        for bw in bws:
            bw.close()

    sys.stderr.write(f"[INFO] Wrote: {args.output}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
