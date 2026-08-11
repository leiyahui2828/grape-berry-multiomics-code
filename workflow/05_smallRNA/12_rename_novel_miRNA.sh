#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

[[ -f scripts/novel_name_mirdp2.pl ]] || {
    echo "[ERROR] Missing scripts/novel_name_mirdp2.pl" >&2
    exit 1
}

for sample in "${SAMPLES[@]}"; do
    perl scripts/novel_name_mirdp2.pl \
      -pbed "novel_miRNA_plant/${sample}.reads/${sample}.reads_filter_P_prediction" \
      -outpre "novel_miRNA_plant/${sample}_novel"
done

perl scripts/novel_name_mirdp2.pl \
  -pbed novel_miRNA_plant/all.reads/all.reads_filter_P_prediction \
  -outpre novel_miRNA_plant/all_novel
