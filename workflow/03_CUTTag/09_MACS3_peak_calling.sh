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

mkdir -p "$OUT/raw_peaks" "$OUT/logs"

# The attached execution log gives exact MACS3-extra parameters only for H3K27ac.
# For other marks, restore the historical config and provide:
#
#   export MACS3_FORMAT=BAM
#   export MACS3_EXTRA='...'
#
if [[ "$MARK" == "H3K27ac" ]]; then
    MACS3_FORMAT="${MACS3_FORMAT:-$H3K27AC_MACS3_FORMAT}"
    MACS3_EXTRA="${MACS3_EXTRA:-$H3K27AC_MACS3_EXTRA}"
else
    if [[ -z "${MACS3_EXTRA:-}" ]]; then
        echo "[ERROR] Exact historical MACS3_EXTRA for $MARK is not present in the uploaded files." >&2
        echo "[INFO] Peak type is already set to: $PEAK_TYPE" >&2
        echo "[INFO] Restore/check the original $MARK config, then run for example:" >&2
        echo "       export MACS3_FORMAT=BAM" >&2
        echo "       export MACS3_EXTRA='--nomodel --shift -100 --extsize 200 --SPMR -B'" >&2
        echo "       bash 09_MACS3_peak_calling.sh $MARK" >&2
        exit 1
    fi
    MACS3_FORMAT="${MACS3_FORMAT:-BAM}"
fi

read -r -a EXTRA_ARR <<< "$MACS3_EXTRA"

BROAD_ARGS=()
if [[ "$PEAK_TYPE" == "broad" ]]; then
    BROAD_ARGS=(--broad)
fi

PLOT_HELPER="scripts/plot_macs_cutoff_analysis.py"

for sample in "${SAMPLES[@]}"; do
    echo "[MACS3:$PEAK_TYPE] $sample"

    macs3 callpeak \
        -g "$EFFECTIVE_GENOME_SIZE" \
        -t "$OUT/clean_bams/${sample}_final.bam" \
        --name "$OUT/raw_peaks/${sample}" \
        --format "$MACS3_FORMAT" \
        --keep-dup all \
        --qvalue 0.075 \
        --cutoff-analysis \
        "${BROAD_ARGS[@]}" \
        "${EXTRA_ARR[@]}" \
        > "$OUT/logs/${sample}_callpeak_with_macs3.log" 2>&1

    if [[ -f "$PLOT_HELPER" ]]; then
        python "$PLOT_HELPER" \
            "$OUT/raw_peaks/${sample}_cutoff_analysis.txt" \
            "$OUT/raw_peaks/${sample}_cutoff_analysis.pdf" \
            > "$OUT/logs/${sample}_plot_macs_cutoff_analysis.log" 2>&1
    fi
done
