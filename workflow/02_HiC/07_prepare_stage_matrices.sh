#!/usr/bin/env bash
set -euo pipefail

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.txt}"

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "[ERROR] Sample list not found: $SAMPLES_FILE" >&2
    exit 1
fi

mapfile -t SAMPLES < <(grep -v '^[[:space:]]*$' "$SAMPLES_FILE")

mkdir -p Matrix/Merged logs

# The supplied analysis contains one Hi-C library per developmental stage.
# Therefore hicSumMatrices reduced to a file copy in the executed dry-run.

for sample in "${SAMPLES[@]}"; do
    stage="${sample%%_Rep*}"

    echo "[Stage matrix] ${sample} -> ${stage}"

    cp \
        "Matrix/Raw/${sample}.cool" \
        "Matrix/Merged/${stage}.cool"

    echo "Copied Matrix/Raw/${sample}.cool to Matrix/Merged/${stage}.cool" \
        > "logs/${stage}_hicSumMatrices.log"
done
