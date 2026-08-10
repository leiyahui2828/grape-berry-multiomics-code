#!/usr/bin/env bash
set -euo pipefail

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.txt}"

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "[ERROR] Sample list not found: $SAMPLES_FILE" >&2
    exit 1
fi

mapfile -t SAMPLES < <(grep -v '^[[:space:]]*$' "$SAMPLES_FILE")

mkdir -p Matrix/Normalized logs

STAGES=()
for sample in "${SAMPLES[@]}"; do
    stage="${sample%%_Rep*}"
    STAGES+=("$stage")
done

INPUTS=()
OUTPUTS=()

for stage in "${STAGES[@]}"; do
    INPUTS+=("Matrix/Merged/${stage}.cool")
    OUTPUTS+=("Matrix/Normalized/${stage}.cool")
done

echo "[hicNormalize] ${STAGES[*]}"

hicNormalize \
    --matrices "${INPUTS[@]}" \
    --outFileName "${OUTPUTS[@]}" \
    --normalize smallest \
    > logs/hicNormalize.log 2>&1
