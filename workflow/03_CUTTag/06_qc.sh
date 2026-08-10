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

mkdir -p "$OUT/quality_control" "$OUT/logs"

COMBINE_FASTP="scripts/combine_fastp_report.R"

if [[ -f "$COMBINE_FASTP" ]]; then
    JSONS=()
    for sample in "${SAMPLES[@]}"; do
        JSONS+=("$OUT/clean_data/${sample}_fastp.json")
    done

    Rscript "$COMBINE_FASTP" \
        --input "${JSONS[@]}" \
        --output "$OUT/quality_control/" \
        > "$OUT/logs/combine_fastp_reports.log" 2>&1
else
    echo "[WARNING] Missing $COMBINE_FASTP; fastp summary step skipped." >&2
fi

BAMS=()
LABELS=()
for sample in "${SAMPLES[@]}"; do
    BAMS+=("$OUT/raw_bams/${sample}_sorted.bam")
    LABELS+=("$sample")
done

BLACKLIST_ARGS=()
if [[ -n "$BLACKLIST" ]]; then
    BLACKLIST_ARGS=(--blackListFileName "$BLACKLIST")
fi

estimateReadFiltering \
    --numberOfProcessors "$FILTER_THREADS" \
    --bam "${BAMS[@]}" \
    --sampleLabels "${LABELS[@]}" \
    --outFile "$OUT/quality_control/estimate_read_filtering.tsv" \
    --ignoreDuplicates \
    --minMappingQuality "$MIN_MAPQ" \
    --samFlagExclude 260 \
    "${BLACKLIST_ARGS[@]}" \
    > "$OUT/logs/estimate_read_filtering.log" 2>&1
