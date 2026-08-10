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

# The H3K27ac source run used IDR.
#
# Other historical mark configs were not supplied. Choose explicitly if needed:
#   export PEAK_SELECTION=idr
# or
#   export PEAK_SELECTION=intersect
#
PEAK_SELECTION="${PEAK_SELECTION:-idr}"
MIN_OVERLAP="${MIN_OVERLAP:-0.5}"

mkdir -p "$OUT/clean_peaks/$PEAK_SELECTION" "$OUT/logs"

for stage in $(get_stages); do
    sample_prefix="${MARK}_${stage}"
    rep1="$OUT/raw_peaks/${sample_prefix}_rep1_peaks.${PEAK_TYPE}Peak"
    rep2="$OUT/raw_peaks/${sample_prefix}_rep2_peaks.${PEAK_TYPE}Peak"
    final="$OUT/clean_peaks/${PEAK_SELECTION}/${sample_prefix}_peaks.${PEAK_TYPE}Peak"

    if [[ "$PEAK_SELECTION" == "idr" ]]; then
        idr_txt="$OUT/clean_peaks/idr/${sample_prefix}_true_rep_idr.txt"

        idr \
            --samples "$rep1" "$rep2" \
            --idr-threshold "$IDR_THRESHOLD" \
            --output-file "$idr_txt" \
            --plot \
            --input-file-type "${PEAK_TYPE}Peak" \
            --rank p.value \
            > "$OUT/logs/${sample_prefix}_select_peaks_by_idr.log" 2>&1

        if [[ "$PEAK_TYPE" == "narrow" ]]; then
            awk -v OFS="\t" 'BEGIN {FS=OFS} {$4=$1":"$2"_"$3; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' \
                "$idr_txt" | sort -k1,1 -k2,2n -k3,3n > "$final"
        else
            # broadPeak contains 9 standard columns.
            awk -v OFS="\t" 'BEGIN {FS=OFS} {$4=$1":"$2"_"$3; print $1,$2,$3,$4,$5,$6,$7,$8,$9}' \
                "$idr_txt" | sort -k1,1 -k2,2n -k3,3n > "$final"
        fi

    elif [[ "$PEAK_SELECTION" == "intersect" ]]; then
        mkdir -p "$OUT/clean_peaks/intersect"

        # The generic source workflow applies reciprocal overlap to q-score-
        # filtered replicate peaks for the intersect route.
        rep1_cut="$OUT/clean_peaks/cutoff/${sample_prefix}_rep1_peaks.${PEAK_TYPE}Peak"
        rep2_cut="$OUT/clean_peaks/cutoff/${sample_prefix}_rep2_peaks.${PEAK_TYPE}Peak"

        bedtools intersect \
            -f "$MIN_OVERLAP" \
            -r \
            -a "$rep1_cut" \
            -b "$rep2_cut" \
            -wa \
            | sort -k1,1 -k2,2n -k3,3n \
            > "$final"
    else
        echo "[ERROR] PEAK_SELECTION must be 'idr' or 'intersect'." >&2
        exit 1
    fi
done

echo "[DONE] Reproducible stage-level peaks: $PEAK_SELECTION"
