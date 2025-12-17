#!/bin/bash

# prepare_cpg_matrix.sh
# Create per-context .tsv.gz matrix files from binom.data_*.txt.gz inputs

usage() {
  echo "Usage: $0 -f <fasta_or_fai> -o <output_dir> --g1 <rep1> [--g1 <rep2> ...] --g2 <rep1> [--g2 <rep2> ...] -W <window_size> -S <step_size> [--verbose]"
  echo ""
  echo "Options:"
  echo "  -f         Reference FASTA or .fai file"
  echo "  -o         Output directory"
  echo "  --g1       Group1 sample directory (can be specified multiple times)"
  echo "  --g2       Group2 sample directory (can be specified multiple times)"
  echo "  -W         Window size (default: 1000)"
  echo "  -S         Step size (default: 500)"
  echo "  --verbose  Print progress information"
  echo "  -h, --help Show this help message"
  exit 1
}

# defaults
WINDOW=1000
STEP=500
VERBOSE=0
THREADS=1

# parse arguments
TEMP=$(getopt -o f:o:W:S:h -l g1:,g2:,verbose,help -n "$0" -- "$@") || usage
eval set -- "$TEMP"

# initialize arrays
G1=()
G2=()

while true; do
  case "$1" in
    -f) REF="$2"; shift 2 ;;
    -o) OUTDIR="$2"; shift 2 ;;
    --g1) G1+=("$2"); shift 2 ;;
    --g2) G2+=("$2"); shift 2 ;;
    -W) WINDOW="$2"; shift 2 ;;
    -S) STEP="$2"; shift 2 ;;
    --verbose) VERBOSE=1; shift ;;
    -h|--help) usage ;;
    --) shift; break ;;
    *) usage ;;
  esac
done

[[ -z "$REF" || -z "$OUTDIR" || ${#G1[@]} -eq 0 || ${#G2[@]} -eq 0 ]] && usage

mkdir -p "$OUTDIR/sample1" "$OUTDIR/sample2" "$OUTDIR/tmp_beds"

# get chromosome length
if [[ "$REF" =~ \.fa(sta)?$ ]]; then
  [[ $VERBOSE -eq 1 ]] && echo "[INFO] Extracting chromosome lengths from $REF..."
  samtools faidx "$REF"
  cut -f1,2 "${REF}.fai" > "$OUTDIR/chrom_lengths.txt"
else
  [[ $VERBOSE -eq 1 ]] && echo "[INFO] Using chromosome lengths file: $REF"
  cp "$REF" "$OUTDIR/chrom_lengths.txt"
fi

# make sliding windows
[[ $VERBOSE -eq 1 ]] && echo "[INFO] Generating sliding windows..."
bedtools makewindows -g "$OUTDIR/chrom_lengths.txt" -w "$WINDOW" -s "$STEP" > "$OUTDIR/tmp_beds/windows.bed"

# function to process one directory
process_sample() {
  local SAMPLE_DIR=$1
  local GROUP=$2
  local OUTDIR=$3

  local SAMPLE_NAME=$(basename "$SAMPLE_DIR")
  for CONTEXT in CpG CHG CHH; do
    [[ $VERBOSE -eq 1 ]] && echo "[INFO] Processing $SAMPLE_DIR for $CONTEXT..."

    gunzip -c "$SAMPLE_DIR"/binom.data_*.txt.gz | \
      awk -v ctx="$CONTEXT" -v group="$GROUP" -v sample="$SAMPLE_NAME" \
      'BEGIN{OFS="\t"} $1 == ctx {print $2,$3,$4,ctx,group,sample,$5,$6,$7,$8}' > "$OUTDIR/tmp_beds/${SAMPLE_NAME}_${CONTEXT}.bed"

    bedtools intersect -a "$OUTDIR/tmp_beds/windows.bed" \
      -b "$OUTDIR/tmp_beds/${SAMPLE_NAME}_${CONTEXT}.bed" -wa -wb | \
      gzip -c > "$OUTDIR/$GROUP/${SAMPLE_NAME}_${CONTEXT}.tsv.gz"

    [[ $VERBOSE -eq 1 ]] && echo "[INFO] Saved $OUTDIR/$GROUP/${SAMPLE_NAME}_${CONTEXT}.tsv.gz"
  done
}

# process group1
for dir in "${G1[@]}"; do
  process_sample "$dir" "sample1" "$OUTDIR"
done

# process group2
for dir in "${G2[@]}"; do
  process_sample "$dir" "sample2" "$OUTDIR"
done

[[ $VERBOSE -eq 1 ]] && echo "[INFO] All processing complete."
