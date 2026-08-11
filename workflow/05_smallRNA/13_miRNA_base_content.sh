#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

[[ -f scripts/mirna_base_content.pl ]] || {
    echo "[ERROR] Missing scripts/mirna_base_content.pl" >&2
    exit 1
}

mkdir -p miRNA_base miRNA_base/novel_bowtie_index logs

for sample in "${SAMPLES[@]}"; do
    awk '$NF=="miRNA"' "reads_mapping/${sample}.out.read_anno.txt" \
      | awk '{print ">"$1"\n"$2}' \
      > "miRNA_base/${sample}.known_mirna.fa"

    perl scripts/mirna_base_content.pl \
      "miRNA_base/${sample}.known_mirna.fa" \
      "miRNA_base/${sample}.known"

    awk '$NF=="unknown" || $NF=="intron_sense" || $NF=="intron_antisense"' \
      "reads_mapping/${sample}.out.read_anno.txt" \
      | awk '{print ">"$1"\n"$2}' \
      > "miRNA_base/${sample}.reads_for_novel.fa"

    cp "novel_miRNA_plant/${sample}_novel.precursor.fa" \
      "miRNA_base/novel_bowtie_index/${sample}_novel.precursor.fa"

    bowtie-build \
      "miRNA_base/novel_bowtie_index/${sample}_novel.precursor.fa" \
      "miRNA_base/novel_bowtie_index/${sample}_novel.precursor"

    bowtie \
      -f -v 0 -p 3 -a --best --strata \
      --al "miRNA_base/${sample}.novel_miRNA.fa" \
      "miRNA_base/novel_bowtie_index/${sample}_novel.precursor" \
      "miRNA_base/${sample}.reads_for_novel.fa" \
      > "miRNA_base/${sample}.novel_miRNA.bwt" \
      2> "logs/${sample}.novel_miRNA_bwt.log"

    perl scripts/mirna_base_content.pl \
      "miRNA_base/${sample}.novel_miRNA.fa" \
      "miRNA_base/${sample}.novel_miRNA"
done

cat miRNA_base/*known_mirna.fa > miRNA_base/all.known_mirna.fa
perl scripts/mirna_base_content.pl \
  miRNA_base/all.known_mirna.fa miRNA_base/all.known_mirna

cat miRNA_base/*novel_miRNA.fa > miRNA_base/all.novel_miRNA.fa
perl scripts/mirna_base_content.pl \
  miRNA_base/all.novel_miRNA.fa miRNA_base/all.novel_miRNA
