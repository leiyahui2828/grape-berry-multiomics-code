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
mkdir -p "$OUT/clean_peaks/merge"

PEAKS=()
for stage in $(get_stages); do
    PEAKS+=("$OUT/clean_peaks/${PEAK_SELECTION}/${MARK}_${stage}_peaks.${PEAK_TYPE}Peak")
done

cat "${PEAKS[@]}" > "$OUT/clean_peaks/merge/cat_peaks.bed"

bedtools sort \
    -i "$OUT/clean_peaks/merge/cat_peaks.bed" \
    -g "${GENOME}.fai" \
    > "$OUT/clean_peaks/merge/cat_sorted_peaks.bed"

bedtools merge \
    -i "$OUT/clean_peaks/merge/cat_sorted_peaks.bed" \
    > "$OUT/clean_peaks/merge/merged_peaks.bed"

awk 'OFS="\t" {print $1":"$2"-"$3,$1,$2+1,$3,"."}' \
    "$OUT/clean_peaks/merge/merged_peaks.bed" \
    > "$OUT/clean_peaks/merge/merged_peaks.saf"

mkdir -p "$OUT/peak_annotation" "$OUT/logs"

uropa \
    --bed "$OUT/clean_peaks/merge/merged_peaks.bed" \
    --gtf "$GTF" \
    --feature gene \
    --feature-anchor start center end \
    --distance 2000 500 \
    --relative-location \
        Upstream Downstream OverlapStart OverlapEnd \
        PeakInsideFeature FeatureInsidePeak \
    --show-attributes gene_id \
    --outdir "$OUT/peak_annotation" \
    > "$OUT/logs/merged_peak_annotation_with_uropa.log" 2>&1
