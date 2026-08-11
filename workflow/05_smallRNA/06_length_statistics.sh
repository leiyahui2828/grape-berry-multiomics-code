#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

for helper in scripts/stat_srna.pl scripts/common_specific_reads.pl; do
    [[ -f "$helper" ]] || { echo "[ERROR] Missing $helper" >&2; exit 1; }
done

for sample in "${SAMPLES[@]}"; do
    perl scripts/stat_srna.pl \
      "clean_data/uniq_data/${sample}.fa" \
      "clean_data/uniq_data/${sample}" \
      > "logs/stat_srna_${sample}.log" 2>&1
done

ls clean_data/uniq_data/*.fa \
  | awk -F "/" '{print $3"\t"$0}' \
  | sed 's/.fa//' \
  > clean_data/uniq_data/clean_srna.list

perl scripts/common_specific_reads.pl \
  clean_data/uniq_data/clean_srna.list \
  > clean_data/uniq_data/clean_srna.count.tsv

cat clean_data/uniq_data/SR_PN*.len.total.txt \
  > clean_data/uniq_data/all.len.total.txt

cat clean_data/uniq_data/SR_PN*.len.uniq.txt \
  > clean_data/uniq_data/all.len.uniq.txt
