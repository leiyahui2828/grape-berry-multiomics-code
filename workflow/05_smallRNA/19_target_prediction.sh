#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh

mkdir -p target_gene_plant logs

cp miRNA_exp/all_know_novel.mature.fa \
  target_gene_plant/all.mature.fa

[[ -f ref/exon.fa ]] || {
    echo "[ERROR] Missing ref/exon.fa" >&2
    exit 1
}
cp ref/exon.fa target_gene_plant/exon.fa

[[ -x "$PSROBOT" || -f "$PSROBOT" ]] || {
    echo "[ERROR] Missing psRobot executable: $PSROBOT" >&2
    exit 1
}

"$PSROBOT" \
  -s target_gene_plant/all.mature.fa \
  -t target_gene_plant/exon.fa \
  -o target_gene_plant/psRobot.out \
  > logs/psRobot_target_gene_predict.log 2>&1

grep "^>" target_gene_plant/psRobot.out \
  | sed 's/>//' \
  | awk 'NF>2 {print $4"\t"$1}' \
  | sort -u \
  > target_gene_plant/psRobot.pair

[[ -f "$TARGETFINDER" || -x "$TARGETFINDER" ]] || {
    echo "[ERROR] Missing TargetFinder wrapper: $TARGETFINDER" >&2
    exit 1
}

"$TARGETFINDER" \
  -f target_gene_plant/all.mature.fa \
  -d target_gene_plant/exon.fa \
  -o target_gene_plant/targetfinder.out \
  -p table \
  > logs/targetfinder_target_gene_predict.log 2>&1

grep -v "^No results" target_gene_plant/targetfinder.out \
  | awk '{print $2"\t"$1}' \
  | sort -u \
  > target_gene_plant/targetfinder.pair

echo "[INFO] psRNATarget was run online in the historical workflow."
echo "[INFO] Place downloaded output at target_gene_plant/psRNATarget.out."
