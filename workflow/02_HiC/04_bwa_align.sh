#!/usr/bin/env bash
set -euo pipefail

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.txt}"

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "[ERROR] Sample list not found: $SAMPLES_FILE" >&2
    exit 1
fi

mapfile -t SAMPLES < <(grep -v '^[[:space:]]*$' "$SAMPLES_FILE")

mkdir -p aligned logs

for sample in "${SAMPLES[@]}"; do
    for mate in R1 R2; do
        echo "[BWA] ${sample} ${mate}"

        bwa mem \
            -A1 \
            -B4 \
            -E50 \
            -L0 \
            -t 20 \
            -o "aligned/${sample}_${mate}.sam" \
            ref/genome \
            "clean_data/${sample}_${mate}.fastq.gz" \
            2> "logs/${sample}_${mate}_bwa.log"
    done
done
