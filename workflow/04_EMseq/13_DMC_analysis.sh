#!/usr/bin/env bash
set -euo pipefail
HELPER="software/run_DM_analysis.R"

for f in "$HELPER" samples_info.txt contrasts.txt; do
    [[ -f "$f" ]] || { echo "[ERROR] Missing $f" >&2; exit 1; }
done

mkdir -p DM_Analysis logs

for context in CpG CHG CHH; do
    Rscript "$HELPER" \
      -s samples_info.txt \
      -c contrasts.txt \
      -f calmeth \
      --context "$context" \
      -t 8 \
      -o "DM_Analysis/${context}_site" \
      --dm dmc \
      > "logs/${context}_site_DM.log" 2>&1
done
