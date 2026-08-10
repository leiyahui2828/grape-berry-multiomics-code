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

mkdir -p "$OUT/clean_bams" "$OUT/logs"

BLACKLIST_ARGS=()
if [[ -n "$BLACKLIST" ]]; then
    BLACKLIST_ARGS=(--blackListFileName "$BLACKLIST")
fi

for sample in "${SAMPLES[@]}"; do
    echo "[alignmentSieve] $sample"
    alignmentSieve \
        --numberOfProcessors "$FILTER_THREADS" \
        --bam "$OUT/raw_bams/${sample}_sorted.bam" \
        --outFile "$OUT/clean_bams/${sample}_final.bam" \
        --filterMetrics "$OUT/logs/${sample}_sieve_alignment.metrics.tsv" \
        --ignoreDuplicates \
        --minMappingQuality "$MIN_MAPQ" \
        --samFlagExclude 260 \
        "${BLACKLIST_ARGS[@]}" \
        > "$OUT/logs/${sample}_sieve_alignment.log" 2>&1

    samtools index "$OUT/clean_bams/${sample}_final.bam"
done
