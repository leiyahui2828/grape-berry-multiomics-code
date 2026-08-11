#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

[[ -f "$MAPPER_PL" ]] || { echo "[ERROR] Missing $MAPPER_PL" >&2; exit 1; }

mkdir -p clean_data/uniq_data logs

i=1
for sample in "${SAMPLES[@]}"; do
    stage="$(get_stage "$sample")"
    rep="${sample##*rep}"

    case "$stage" in
        PN4)  prefix="P1${rep}" ;;
        PN7)  prefix="P2${rep}" ;;
        PN10) prefix="P3${rep}" ;;
        PN15) prefix="P4${rep}" ;;
        PN18) prefix="P5${rep}" ;;
        *) echo "[ERROR] Unknown stage: $stage" >&2; exit 1 ;;
    esac

    perl "$MAPPER_PL" \
      "clean_data/${sample}_R1_trimmed.fq" \
      -e \
      -g "$prefix" \
      -h \
      -m \
      -s "clean_data/uniq_data/${sample}.fa" \
      > "logs/${sample}_fqtofa.log" 2>&1
done
