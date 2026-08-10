#!/usr/bin/env bash
set -euo pipefail

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.txt}"

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "[ERROR] Sample list not found: $SAMPLES_FILE" >&2
    exit 1
fi

mapfile -t SAMPLES < <(grep -v '^[[:space:]]*$' "$SAMPLES_FILE")

mkdir -p Matrix/Final logs

RESOLUTIONS="1000,2000,5000,10000,15000,20000,25000,40000,100000,250000,500000,1000000"

STAGES=()
for sample in "${SAMPLES[@]}"; do
    stage="${sample%%_Rep*}"
    STAGES+=("$stage")
done

for stage in "${STAGES[@]}"; do
    echo "[cooler zoomify] ${stage}"

    cooler zoomify \
        "Matrix/Normalized/${stage}.cool" \
        -n 4 \
        --resolutions "$RESOLUTIONS" \
        -o "Matrix/Final/${stage}_nobalanced.mcool" \
        > "logs/${stage}_cooler_zoomify_nobalanced.log" 2>&1

    cooler zoomify \
        "Matrix/Normalized/${stage}.cool" \
        -n 4 \
        --balance \
        --resolutions "$RESOLUTIONS" \
        -o "Matrix/Final/${stage}.mcool" \
        > "logs/${stage}_cooler_zoomify_balanced.log" 2>&1
done
