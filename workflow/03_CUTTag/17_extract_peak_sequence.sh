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

PEAK_SELECTION="${PEAK_SELECTION:-idr}"
TOP_N="${TOP_N:-5000}"
SUMMIT_FLANK="${SUMMIT_FLANK:-100}"

mkdir -p "$OUT/clean_peaks/${PEAK_SELECTION}"

for stage in $(get_stages); do
    sample_prefix="${MARK}_${stage}"
    peak="$OUT/clean_peaks/${PEAK_SELECTION}/${sample_prefix}_peaks.${PEAK_TYPE}Peak"
    peak_fa="$OUT/clean_peaks/${PEAK_SELECTION}/${sample_prefix}_peaks.fa"

    sort -k8,8nr "$peak" \
        | head -n "$TOP_N" \
        | bedtools getfasta -fi "$GENOME" -bed - \
        > "$peak_fa"

    if [[ "$PEAK_TYPE" == "narrow" ]]; then
        summit_fa="$OUT/clean_peaks/${PEAK_SELECTION}/${sample_prefix}_summit.fa"

        sort -k8,8nr "$peak" \
            | head -n "$TOP_N" \
            | awk -v OFS='\t' -v flank="$SUMMIT_FLANK" \
                '{print $1,$2+$10-flank,$2+$10+flank+1}' \
            | awk '$2>=0' \
            | bedtools getfasta -fi "$GENOME" -bed - \
            > "$summit_fa"
    fi
done
