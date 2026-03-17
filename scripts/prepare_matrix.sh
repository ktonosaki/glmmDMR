#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

show_help() {
    echo "Usage: prepare_matrix_group_based.sh --fasta <genome.fa or .fai> \\"
    echo "                             --group1 <file1.tsv.gz> <file2.tsv.gz> ... \\"
    echo "                             --group2 <file1.tsv.gz> <file2.tsv.gz> ... \\"
    echo "                             [--group_labels group1_label group2_label] \\"
    echo "                             [--window 1000] \\"
    echo "                             [--slide 500] \\"
    echo "                             [--output outdir] \\"
    echo "                             [--tmpdir tmpdir] \\"
    echo "                             [--help]"
    echo ""
    echo "Options:"
    echo "  --fasta            Fasta or .fai file used to generate genome windows"
    echo "  --group1           List of binomtest_result.tsv.gz files for group 1"
    echo "  --group2           List of binomtest_result.tsv.gz files for group 2"
    echo "  --group_labels     Labels for group1 and group2 (default: group1 group2)"
    echo "  --window           Sliding window size in bp (default: 1000)"
    echo "  --slide            Sliding window step in bp (default: 500)"
    echo "  --output           Output directory (default: ./matrix_out)"
    echo "  --tmpdir           Temporary directory for intermediate files"
    echo "  --help             Show this message and exit"
    exit 0
}

# Default values
window=300
slide=200
outdir="./matrix_out"
tmpdir=""
group1_label="group1"
group2_label="group2"
group1_files=()
group2_files=()

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta) fasta="$2"; shift 2 ;;
        --group1) shift; while [[ $# -gt 0 && "$1" != --* ]]; do group1_files+=("$1"); shift; done ;;
        --group2) shift; while [[ $# -gt 0 && "$1" != --* ]]; do group2_files+=("$1"); shift; done ;;
        --group_labels) group1_label="$2"; group2_label="$3"; shift 3 ;;
        --window) window="$2"; shift 2 ;;
        --slide) slide="$2"; shift 2 ;;
        --output) outdir="$2"; shift 2 ;;
        --tmpdir) tmpdir="$2"; shift 2 ;;
        --help) show_help ;;
        *) echo "Unknown option $1"; exit 1 ;;
    esac
done

# Check required options
if [[ -z "$fasta" || ${#group1_files[@]} -eq 0 || ${#group2_files[@]} -eq 0 ]]; then
    echo "Error: --fasta, --group1, and --group2 are required."
    echo "Use --help to see usage."
    exit 1
fi

# Validate input files
for f in "${group1_files[@]}" "${group2_files[@]}"; do
    [[ -f "$f" ]] || { echo "Error: input file not found: $f"; exit 1; }
done

# Compute common prefix/suffix across all sample basenames
all_bases=()
for f in "${group1_files[@]}" "${group2_files[@]}"; do
    all_bases+=("$(basename "$f" .tsv.gz)")
done

common_prefix() {
    local arr=("$@")
    local prefix="${arr[0]}"
    for s in "${arr[@]}"; do
        while [[ -n "$prefix" && "${s#"$prefix"}" == "$s" ]]; do
            prefix="${prefix%?}"
        done
    done
    echo "$prefix"
}

common_suffix() {
    local arr=("$@")
    local suffix="${arr[0]}"
    for s in "${arr[@]}"; do
        while [[ -n "$suffix" && "${s%"$suffix"}" == "$s" ]]; do
            suffix="${suffix#?}"
        done
    done
    echo "$suffix"
}

common_pref=$(common_prefix "${all_bases[@]}")
common_suf=$(common_suffix "${all_bases[@]}")
echo "[INFO] Common sample prefix: '${common_pref}'"
echo "[INFO] Common sample suffix: '${common_suf}'"

# Prepare output and temp dirs
mkdir -p "$outdir"
if [[ -z "$tmpdir" ]]; then
    tmpdir=$(mktemp -d)
    trap 'rm -rf "$tmpdir"' EXIT
else
    mkdir -p "$tmpdir"
fi

# Check required tools
command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "Error: samtools not found"; exit 1; }
command -v zcat >/dev/null 2>&1 || { echo "Error: zcat not found"; exit 1; }

# Prepare .fai if needed
if [[ "$fasta" == *.fai ]]; then
    fai="$fasta"
else
    fai="${fasta}.fai"
    [[ -f "$fai" ]] || samtools faidx "$fasta"
fi

# Create window bed file
window_bed="$tmpdir/windows.bed"
bedtools makewindows -g "$fai" -w "$window" -s "$slide" > "$window_bed"
if [[ ! -s "$window_bed" ]]; then
    echo "Error: window bed is empty. Check fasta/fai and window/slide settings."
    exit 1
fi

# Create BED from each input
process_group() {
    label="$1"
    context="$2"
    shift 2
    files=("$@")
    for f in "${files[@]}"; do
        base=$(basename "$f" .tsv.gz)
        sample_id="$base"
        if [[ -n "$common_pref" ]]; then
            sample_id="${sample_id#"$common_pref"}"
        fi
        if [[ -n "$common_suf" ]]; then
            sample_id="${sample_id%"$common_suf"}"
        fi
        if [[ -z "$sample_id" ]]; then
            sample_id="$base"
        fi
        echo "[INFO] Processing $context for $label: $sample_id"
        zcat "$f" | awk -v OFS='\t' -v c="$context" -v g="$label" -v s="$sample_id" '
            NR>1 && $6 == c {
                cov = $4 + $5;
                print $1, $2, $2+1, $6, g, s, $3, $4, $5, cov
            }
        ' > "$tmpdir/${label}_${base}_${context}_cytosines.bed"
    done
}

# Main loop for CpG, CHG, CHH
for ctx in CpG CHG CHH; do
    echo "[INFO] Processing context: $ctx"

    # process both groups
    process_group "$group1_label" "$ctx" "${group1_files[@]}"
    process_group "$group2_label" "$ctx" "${group2_files[@]}"

    mapfile -t all_bed < <(find "$tmpdir" -name "*_${ctx}_cytosines.bed" -print)
    if [[ ${#all_bed[@]} -eq 0 ]]; then
        echo "[WARNING] No $ctx cytosine BED files found"
        continue
    fi

    cat "${all_bed[@]}" | LC_ALL=C sort -k1,1 -k2,2n | bedtools intersect -a "$window_bed" -b - -wa -wb |
    awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' |
    gzip > "$outdir/${group1_label}_${group2_label}_${ctx}_matrix.tsv.gz"

    # Cleanup per-context temp beds to save space
    rm -f "${all_bed[@]}"
done

echo "[DONE] Matrix files written to $outdir"
