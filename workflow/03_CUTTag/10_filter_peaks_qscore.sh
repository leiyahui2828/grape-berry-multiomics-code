#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: bash $0 <H3K27ac|H3K4me3|H3K9ac|H3K27me3>" >&2
    exit 1
fi

MARK="$1"
source ./00_config.sh
PEAK_TYPE="$(get_peak_type "$MARK")"
OUT="$(result_dir "$MARK")"
mapfile -t SAMPLES < <(get_samples "$MARK")

if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
    echo "[ERROR] No samples found for mark: $MARK" >&2
    exit 1
fi

mkdir -p "$OUT/clean_peaks/cutoff"

for sample in "${SAMPLES[@]}"; do
    infile="$OUT/raw_peaks/${sample}_peaks.${PEAK_TYPE}Peak"
    outfile="$OUT/clean_peaks/cutoff/${sample}_peaks.${PEAK_TYPE}Peak"

    echo "[q-score > $QSCORE_CUTOFF] $sample"
    awk -v q="$QSCORE_CUTOFF" '$9>q' "$infile" > "$outfile"
done
