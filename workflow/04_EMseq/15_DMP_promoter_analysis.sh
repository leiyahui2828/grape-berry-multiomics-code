#!/usr/bin/env bash
set -euo pipefail
HELPER="software/run_DM_analysis.R"

for f in "$HELPER" samples_info.txt contrasts.txt ref/genes.gff; do
    [[ -f "$f" ]] || { echo "[ERROR] Missing $f" >&2; exit 1; }
done

mkdir -p DM_Analysis logs

for context in CpG CHG CHH; do
    Rscript "$HELPER" \
      -s samples_info.txt \
      -c contrasts.txt \
      -f calmeth \
      --context "$context" \
      --gff ref/genes.gff \
      -t 5 \
      -o "DM_Analysis/${context}_promoter" \
      --dm dmp \
      > "logs/${context}_pro_DM.log" 2>&1
done
