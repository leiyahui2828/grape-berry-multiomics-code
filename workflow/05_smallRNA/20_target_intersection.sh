#!/usr/bin/env bash
set -euo pipefail

for f in \
    scripts/map_data_ids \
    ref/genome.gff3 \
    target_gene_plant/targetfinder.pair \
    target_gene_plant/psRobot.pair \
    target_gene_plant/psRNATarget.out
do
    [[ -f "$f" ]] || { echo "[ERROR] Missing $f" >&2; exit 1; }
done

awk 'NF>1 {print $2"\t"$1}' \
  target_gene_plant/psRNATarget.out \
  | sed 1d \
  > target_gene_plant/psRNATarget.pair

awk '$3=="mRNA"' ref/genome.gff3 \
  | cut -f 9 \
  | awk -F";" '{print $1"\t"$2}' \
  | sed 's/Parent=// ; s/ID=//' \
  > target_gene_plant/trans2gene.txt

for method in targetfinder psRobot psRNATarget; do
    perl scripts/map_data_ids \
      target_gene_plant/trans2gene.txt \
      "target_gene_plant/${method}.pair" \
      -col 1
done

awk '{print $1}' target_gene_plant/psRNATarget.pair \
  | sort -u > target_gene_plant/psRNATarget.pair.uniq

awk '{print $1}' target_gene_plant/targetfinder.pair \
  | sort -u > target_gene_plant/targetfinder.pair.uniq

awk '{print $1}' target_gene_plant/psRobot.pair \
  | sort -u > target_gene_plant/psRobot.pair.uniq

cat \
  target_gene_plant/targetfinder.pair.uniq \
  target_gene_plant/psRobot.pair.uniq \
  target_gene_plant/psRNATarget.pair.uniq \
  | sort \
  | uniq -c \
  > target_gene_plant/target_mirna.all.out

awk '$1==3 {print $2"\t"$3}' \
  target_gene_plant/target_mirna.all.out \
  > target_gene_plant/target_mirna.out
