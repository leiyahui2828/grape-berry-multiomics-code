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

HELPER="software/RunFeatureCounts/abundance_estimates_to_matrix.pl"

if [[ ! -f "$HELPER" ]]; then
    echo "[ERROR] Missing original helper: $HELPER" >&2
    exit 1
fi

mkdir -p "$OUT/counts" "$OUT/logs"

: > "$OUT/counts/count.fofn"
for sample in "${SAMPLES[@]}"; do
    echo "$OUT/counts/${sample}.count" >> "$OUT/counts/count.fofn"
done

perl "$HELPER" \
    --est_method featureCounts \
    --quant_files "$OUT/counts/count.fofn" \
    --out_prefix "$OUT/counts/merged_peaks" \
    > "$OUT/logs/create_count_matrix.log" 2>&1
