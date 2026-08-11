#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

[[ -f scripts/blastRfam2family.pl ]] || {
    echo "[ERROR] Missing scripts/blastRfam2family.pl" >&2
    exit 1
}

mkdir -p miRNA_family

for sample in "${SAMPLES[@]}"; do
    cat \
      "$MIRBASE_DIR/hairpin.${SPECIES_PREFIX}.fa" \
      "novel_miRNA_plant/${sample}_novel.precursor.fa" \
      > "miRNA_family/${sample}_all.hairpin.fa"

    blastn \
      -db "$RFAM_DIR/Rfam_MIR.fasta" \
      -query "miRNA_family/${sample}_all.hairpin.fa" \
      -out "miRNA_family/${sample}_all.hairpin.Rfam_MIR.blast" \
      -outfmt 6 \
      -evalue 0.01

    perl scripts/blastRfam2family.pl \
      "miRNA_family/${sample}_all.hairpin.Rfam_MIR.blast" \
      "$RFAM_DIR/family.txt.MIR" \
      > "miRNA_family/${sample}_all.hairpin.family.out"
done

cat \
  "$MIRBASE_DIR/hairpin.${SPECIES_PREFIX}.fa" \
  novel_miRNA_plant/all_novel.precursor.fa \
  > miRNA_family/all.hairpin.fa

blastn \
  -db "$RFAM_DIR/Rfam_MIR.fasta" \
  -query miRNA_family/all.hairpin.fa \
  -out miRNA_family/all.hairpin.Rfam_MIR.blast \
  -outfmt 6 \
  -evalue 0.01

perl scripts/blastRfam2family.pl \
  miRNA_family/all.hairpin.Rfam_MIR.blast \
  "$RFAM_DIR/family.txt.MIR" \
  > miRNA_family/all.hairpin.family.out
