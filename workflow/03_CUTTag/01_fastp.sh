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

mkdir -p "$OUT/clean_data" "$OUT/logs"

for sample in "${SAMPLES[@]}"; do
    echo "[fastp] $sample"
    fastp \
        -i "raw_data/${sample}_R1.fastq.gz" \
        -o "$OUT/clean_data/${sample}_R1.fastq.gz" \
        -I "raw_data/${sample}_R2.fastq.gz" \
        -O "$OUT/clean_data/${sample}_R2.fastq.gz" \
        --html "$OUT/clean_data/${sample}_fastp.html" \
        --json "$OUT/clean_data/${sample}_fastp.json" \
        --thread "$FASTP_THREADS" \
        --qualified_quality_phred 20 \
        --length_required 50 \
        > "$OUT/logs/${sample}_fastp_quality_control.log" 2>&1
done
