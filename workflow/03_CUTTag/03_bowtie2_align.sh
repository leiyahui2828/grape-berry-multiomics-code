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

mkdir -p "$OUT/raw_bams" "$OUT/logs"

for sample in "${SAMPLES[@]}"; do
    echo "[Bowtie2] $sample"
    bowtie2 \
        -p "$BOWTIE2_THREADS" \
        -x ref/genome \
        -1 "$OUT/clean_data/${sample}_R1.fastq.gz" \
        -2 "$OUT/clean_data/${sample}_R2.fastq.gz" \
        -S "$OUT/raw_bams/${sample}.sam" \
        --no-mixed \
        --no-discordant \
        > "$OUT/logs/${sample}_bowtie2_align.log" 2>&1
done
