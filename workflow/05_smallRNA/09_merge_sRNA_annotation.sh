#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

[[ -f scripts/srna_anno.pl ]] || { echo "[ERROR] Missing scripts/srna_anno.pl" >&2; exit 1; }

for sample in "${SAMPLES[@]}"; do
    perl scripts/srna_anno.pl \
      -fa "reads_mapping/${sample}.mapped.fa" \
      -rfam "$RFAM_DIR/family.txt1" \
      -mirna "reads_mapping/${sample}.miRNA.bwt" \
      -ncrna "reads_mapping/${sample}.ncRNA.bwt" \
      -intron "reads_mapping/${sample}.intron.bwt" \
      -repeat "reads_mapping/${sample}.repeat.bwt" \
      -exon "reads_mapping/${sample}.exon.bwt" \
      -outpre "reads_mapping/${sample}.out" \
      > "logs/${sample}.merge_bowtie.log" 2>&1
done

echo -e "sample\tann_type\ttotal_reads\ttotal_percent\tuniq_reads\tuniq_percent" \
  > reads_mapping/merged_anno_stat_top6.tsv

for file in reads_mapping/*.out.read_anno.stat; do
    sample="$(basename "$file" | cut -d'.' -f1)"
    awk -v s="$sample" 'NR>1 && NR<=7 {print s"\t"$0}' "$file" \
      >> reads_mapping/merged_anno_stat_top6.tsv
done
