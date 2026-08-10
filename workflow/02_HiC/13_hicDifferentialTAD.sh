#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 13_hicDifferentialTAD.sh
# Pairwise differential TAD analysis using HiCExplorer hicDifferentialTAD
#
# Parameters preserved from the supplied hic.code:
#   resolution = 40000
#   --pValue 0.05
#   --mode all
#   --modeReject one
#   --threads 4
#
# For each pair A vs B, the TAD domains called in A are used as --tadDomains,
# matching the original workflow.
###############################################################################

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.txt}"
RESOLUTION=40000
THREADS=4

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "[ERROR] Sample list not found: $SAMPLES_FILE" >&2
    exit 1
fi

mapfile -t SAMPLES < <(grep -v '^[[:space:]]*$' "$SAMPLES_FILE")

STAGES=()
for sample in "${SAMPLES[@]}"; do
    stage="${sample%%_Rep*}"
    STAGES+=("$stage")
done

mkdir -p "TADs/HiCExplorer/${RESOLUTION}" logs

n=${#STAGES[@]}

for ((i=0; i<n; i++)); do
    for ((j=i+1; j<n; j++)); do
        target="${STAGES[$i]}"
        control="${STAGES[$j]}"

        target_matrix="Matrix/Final/${target}.mcool::/resolutions/${RESOLUTION}"
        control_matrix="Matrix/Final/${control}.mcool::/resolutions/${RESOLUTION}"
        tad_domains="TADs/HiCExplorer/${RESOLUTION}/${target}_domains.bed"
        outprefix="TADs/HiCExplorer/${RESOLUTION}/${target}_vs_${control}"
        logfile="logs/${target}_vs_${control}_${RESOLUTION}_hicDifferentialTAD.log"

        for f in "Matrix/Final/${target}.mcool" \
                 "Matrix/Final/${control}.mcool" \
                 "$tad_domains"; do
            if [[ ! -f "$f" ]]; then
                echo "[ERROR] Missing input: $f" >&2
                exit 1
            fi
        done

        echo "[hicDifferentialTAD] ${target} vs ${control}"

        hicDifferentialTAD \
            --targetMatrix "$target_matrix" \
            --controlMatrix "$control_matrix" \
            --tadDomains "$tad_domains" \
            --outFileNamePrefix "$outprefix" \
            --pValue 0.05 \
            --mode all \
            --modeReject one \
            --threads "$THREADS" \
            > "$logfile" 2>&1
    done
done

echo "[DONE] Pairwise differential TAD analysis completed."
