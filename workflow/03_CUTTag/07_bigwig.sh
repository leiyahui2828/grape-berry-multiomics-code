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

for sample in "${SAMPLES[@]}"; do
    echo "[bamCoverage] $sample"
    bamCoverage \
        -p "$BAMCOVERAGE_THREADS" \
        --bam "$OUT/clean_bams/${sample}_final.bam" \
        -o "$OUT/clean_bams/${sample}.bw" \
        --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE" \
        --binSize 10 \
        --normalizeUsing RPKM \
        > "$OUT/logs/${sample}_convert_bam_to_bigwig.log" 2>&1
done
