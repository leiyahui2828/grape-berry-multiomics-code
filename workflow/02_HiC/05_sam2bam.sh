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
        echo "[SAMtools view] ${sample} ${mate}"

        samtools view \
            -hb \
            --threads 6 \
            -o "aligned/${sample}_${mate}.bam" \
            "aligned/${sample}_${mate}.sam" \
            > "logs/${sample}_${mate}_sam2bam.log" 2>&1
    done
done
