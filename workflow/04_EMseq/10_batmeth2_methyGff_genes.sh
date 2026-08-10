#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

mkdir -p FeatureML logs
for sample in "${SAMPLES[@]}"; do
    run_meth_tools BatMeth2 methyGff \
      --coverage "$METH_COVERAGE" \
      -nC "$METH_NC" \
      --genome "$GENOME" \
      -b5 ref/genes.bed5 \
      --body --promoter --GENE \
      --distance "$FEATURE_DISTANCE" \
      --step "$FEATURE_STEP" \
      --methratio "site_samll_window_ML/${sample}.methratio.txt" \
      --out "FeatureML/${sample}.genes" \
      > "logs/${sample}.BatMeth2.gene.log" 2>&1
done
