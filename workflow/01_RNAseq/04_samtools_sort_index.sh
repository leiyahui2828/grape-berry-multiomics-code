#!/usr/bin/env bash
set -euo pipefail

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.txt}"

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "[ERROR] Sample list not found: $SAMPLES_FILE" >&2
    exit 1
fi

mapfile -t SAMPLES < <(grep -v '^\s*$' "$SAMPLES_FILE")

mkdir -p aligned logs

for sample in "${SAMPLES[@]}"; do
    echo "[SAMtools] ${sample}"
    samtools sort \
        -@ 4 \
        -o "aligned/${sample}.bam" \
        "aligned/${sample}.sam" \
        > "logs/${sample}.sort.log" 2>&1

    samtools index "aligned/${sample}.bam"
done
