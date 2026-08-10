#!/usr/bin/env bash
set -euo pipefail

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.txt}"

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "[ERROR] Sample list not found: $SAMPLES_FILE" >&2
    exit 1
fi

mapfile -t SAMPLES < <(grep -v '^\s*$' "$SAMPLES_FILE")

mkdir -p clean_data logs

for sample in "${SAMPLES[@]}"; do
    echo "[fastp] ${sample}"
    fastp \
        -i "raw_data/${sample}_1.fastq.gz" \
        -I "raw_data/${sample}_2.fastq.gz" \
        -o "clean_data/${sample}_1.fastp.fastq.gz" \
        -O "clean_data/${sample}_2.fastp.fastq.gz" \
        -h "clean_data/${sample}.fastp.html" \
        -j "clean_data/${sample}.fastp.json" \
        -w 4 \
        > "logs/${sample}.fastp.log" 2>&1
done
