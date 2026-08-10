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

mkdir -p "$OUT/peak_annotation" "$OUT/logs"

for stage in $(get_stages); do
    sample_prefix="${MARK}_${stage}"
    peak="$OUT/clean_peaks/${PEAK_SELECTION}/${sample_prefix}_peaks.${PEAK_TYPE}Peak"

    echo "[UROPA] $sample_prefix"

    uropa \
        --bed "$peak" \
        --gtf "$GTF" \
        --feature gene \
        --feature-anchor start center end \
        --distance 2000 500 \
        --relative-location \
            Upstream Downstream OverlapStart OverlapEnd \
            PeakInsideFeature FeatureInsidePeak \
        --show-attributes gene_id \
        --outdir "$OUT/peak_annotation" \
        > "$OUT/logs/${sample_prefix}_peak_annotation_with_uropa.log" 2>&1
done
