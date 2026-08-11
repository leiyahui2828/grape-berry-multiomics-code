#!/usr/bin/env bash
set -euo pipefail

for f in \
    scripts/run_DE_analysis.pl \
    miRNA_exp/miRNAs_expressed.count.txt \
    miRNA_exp/sample.txt \
    miRNA_exp/contrasts.txt
do
    [[ -f "$f" ]] || { echo "[ERROR] Missing $f" >&2; exit 1; }
done

mkdir -p DE_result logs

perl scripts/run_DE_analysis.pl \
  --matrix miRNA_exp/miRNAs_expressed.count.txt \
  --samples_file miRNA_exp/sample.txt \
  --contrasts miRNA_exp/contrasts.txt \
  --output DE_result \
  --method DESeq2 \
  > logs/DE_analysis.log 2>&1

if [[ -f scripts/process_files.sh ]]; then
    bash scripts/process_files.sh
else
    echo "[WARNING] scripts/process_files.sh missing; modification step skipped." >&2
fi

awk 'FNR==1 && NR!=1 {next} {print}' \
  DE_result/*DESeq2.modified.DE_results \
  > DE_result/SR_de_result.tsv
