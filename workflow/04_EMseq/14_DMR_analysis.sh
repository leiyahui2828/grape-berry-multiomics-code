#!/usr/bin/env bash
set -euo pipefail
HELPER="software/run_DM_analysis.R"

for f in "$HELPER" samples_info.txt contrasts.txt; do
    [[ -f "$f" ]] || { echo "[ERROR] Missing $f" >&2; exit 1; }
done

mkdir -p DM_Analysis logs

Rscript "$HELPER" -s samples_info.txt -c contrasts.txt -f calmeth \
  --context CpG --window 600 --step 200 -t 5 \
  -o DM_Analysis/CpG_region --dm dmr > logs/CpG_region_DM.log 2>&1

Rscript "$HELPER" -s samples_info.txt -c contrasts.txt -f calmeth \
  --context CHG --window 600 --step 200 -t 8 \
  -o DM_Analysis/CHG_region --dm dmr > logs/CHG_region_DM.log 2>&1

Rscript "$HELPER" -s samples_info.txt -c contrasts.txt -f calmeth \
  --context CHH --window 600 --step 200 -t 5 \
  -o DM_Analysis/CHH_region --dm dmr > logs/CHH_region_DM.log 2>&1
