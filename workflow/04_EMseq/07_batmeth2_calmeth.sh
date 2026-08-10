#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

mkdir -p site_samll_window_ML logs
for sample in "${SAMPLES[@]}"; do
    run_meth_tools BatMeth2 calmeth \
      -Q "$METH_BASEQ" \
      --remove_dup \
      --coverage "$METH_COVERAGE" \
      -nC "$METH_NC" \
      --Regions "$METH_REGION" \
      --step "$METH_STEP" \
      --genome "$GENOME" \
      --binput "raw_bam/${sample}.sorted.bam" \
      --methratio "site_samll_window_ML/${sample}" \
      > "logs/${sample}.BatMeth2.log" 2>&1
done
