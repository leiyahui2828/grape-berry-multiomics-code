#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh

[[ -f "$QUANTIFIER_PL" ]] || { echo "[ERROR] Missing $QUANTIFIER_PL" >&2; exit 1; }

mkdir -p miRNA_exp logs

cat \
  "$MIRBASE_DIR/mature.${SPECIES_PREFIX}.fa" \
  novel_miRNA_plant/all_novel.mature.fa \
  > miRNA_exp/all_know_novel.mature.fa

cat \
  "$MIRBASE_DIR/hairpin.${SPECIES_PREFIX}.fa" \
  novel_miRNA_plant/all_novel.precursor.fa \
  > miRNA_exp/all_know_novel.precursor.fa

cat clean_data/uniq_data/SR_PN*.fa \
  > miRNA_exp/all.reads.fa

cat reads_mapping/SR_PN*.mapped.fa \
  > miRNA_exp/all.mapped.fa

(
  cd miRNA_exp

  perl "../$QUANTIFIER_PL" \
    -p all_know_novel.precursor.fa \
    -m all_know_novel.mature.fa \
    -r all.reads.fa \
    -g 0 \
    > ../logs/all_read_quantifier.log 2>&1

  perl "../$QUANTIFIER_PL" \
    -p all_know_novel.precursor.fa \
    -m all_know_novel.mature.fa \
    -r all.mapped.fa \
    -g 0 \
    > ../logs/all_mapped_read_quantifier.log 2>&1

  # Historical parsing logic.
  cat miRNAs_expressed_all_samples_*.csv \
    | awk '{for(i=5;i<=4+(NF-4)/2;i++){c=c"\t"$i}; print $1 c; c=""}' \
    | awk '{if(NR==1 || !($1 in A)){print $0}; A[$1]=1}' \
    > miRNAs_expressed.count.txt

  cat miRNAs_expressed_all_samples_*.csv \
    | awk '{for(i=5+(NF-4)/2;i<=NF;i++){c=c"\t"$i}; print $1 c; c=""}' \
    | awk '{if(NR==1 || !($1 in A)){print $0}; A[$1]=1}' \
    | sed 's/(norm)//g' \
    > miRNAs_expressed.TPM.txt
)
