#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 12_hicFindTADs.sh
# TAD-like domain calling using HiCExplorer hicFindTADs
#
# Parameters preserved from the supplied hic.code:
#   resolution = 40000
#   --correctForMultipleTesting fdr
#   --thresholdComparisons 0.01
#   --delta 0.01
#   -p 4
###############################################################################

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.txt}"
RESOLUTION=40000
THREADS=4

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "[ERROR] Sample list not found: $SAMPLES_FILE" >&2
    exit 1
fi

mapfile -t SAMPLES < <(grep -v '^[[:space:]]*$' "$SAMPLES_FILE")

mkdir -p "TADs/HiCExplorer/${RESOLUTION}" logs

for sample in "${SAMPLES[@]}"; do
    stage="${sample%%_Rep*}"
    matrix="Matrix/Final/${stage}.mcool::/resolutions/${RESOLUTION}"
    prefix="TADs/HiCExplorer/${RESOLUTION}/${stage}"
    logfile="logs/${stage}_${RESOLUTION}_hicFindTADs.log"

    # Check the mcool container itself; the resolution group is checked by hicFindTADs.
    if [[ ! -f "Matrix/Final/${stage}.mcool" ]]; then
        echo "[ERROR] Missing matrix: Matrix/Final/${stage}.mcool" >&2
        exit 1
    fi

    echo "[hicFindTADs] ${stage}"

    hicFindTADs \
        -p "$THREADS" \
        --matrix "$matrix" \
        --outPrefix "$prefix" \
        --correctForMultipleTesting fdr \
        --thresholdComparisons 0.01 \
        --delta 0.01 \
        > "$logfile" 2>&1
done

echo "[DONE] TAD calling completed."
