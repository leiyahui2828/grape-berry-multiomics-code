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

HELPER="software/RunFeatureCounts/run-featurecounts.R"

if [[ ! -f "$HELPER" ]]; then
    echo "[ERROR] Missing original helper: $HELPER" >&2
    exit 1
fi

mkdir -p "$OUT/counts" "$OUT/logs"

for sample in "${SAMPLES[@]}"; do
    echo "[featureCounts] $sample"
    Rscript "$HELPER" \
        -b "$OUT/clean_bams/${sample}_final.bam" \
        --saf "$OUT/clean_peaks/merge/merged_peaks.saf" \
        --isPairedEnd TRUE \
        -o "$OUT/counts/${sample}" \
        > "$OUT/logs/${sample}_run_feature_counts.log" 2>&1
done
