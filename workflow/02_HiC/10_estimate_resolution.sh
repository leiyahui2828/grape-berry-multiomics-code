#!/usr/bin/env bash
set -euo pipefail

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.txt}"

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "[ERROR] Sample list not found: $SAMPLES_FILE" >&2
    exit 1
fi

mapfile -t SAMPLES < <(grep -v '^[[:space:]]*$' "$SAMPLES_FILE")

HICRES="scripts/hicres.py"

if [[ ! -f "$HICRES" ]]; then
    echo "[ERROR] Missing original helper script: $HICRES" >&2
    echo "[INFO] Copy the exact hicres.py used in the manuscript analysis." >&2
    exit 1
fi

STAGES=()
for sample in "${SAMPLES[@]}"; do
    stage="${sample%%_Rep*}"
    STAGES+=("$stage")
done

for stage in "${STAGES[@]}"; do
    echo "[Hi-C resolution] ${stage}"

    python "$HICRES" \
        "Matrix/Final/${stage}.mcool" \
        --threads 4 \
        --output "Matrix/Final/${stage}_res.csv" \
        --resolutions 1000,2000,5000,10000
done
