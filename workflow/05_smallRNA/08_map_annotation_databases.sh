#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

mkdir -p reads_mapping logs

for sample in "${SAMPLES[@]}"; do
    # Known grape miRNA precursor.
    bowtie \
      -f -v 0 -p 8 -a --best --strata \
      "$MIRBASE_DIR/hairpin.${SPECIES_PREFIX}" \
      "reads_mapping/${sample}.mapped.fa" \
      > "reads_mapping/${sample}.miRNA.bwt" \
      2> "logs/${sample}.miRNA_bowtie.log"

    # Rfam non-miRNA ncRNA.
    bowtie \
      -f -v 1 -p 8 -a --best --strata \
      "$RFAM_DIR/Rfam_rmMIR.fasta" \
      "reads_mapping/${sample}.mapped.fa" \
      > "reads_mapping/${sample}.ncRNA.bwt" \
      2> "logs/${sample}.ncRNA_bowtie.log"

    # Historical commands below use Bowtie1 against indexes named repeat/intron/exon.
    # See METHOD_AUDIT.md for the Bowtie1/Bowtie2 index inconsistency.
    bowtie \
      -f -v 0 -p 2 -a --best --strata \
      ref/repeat \
      "reads_mapping/${sample}.mapped.fa" \
      > "reads_mapping/${sample}.repeat.bwt" \
      2> "logs/${sample}.repeat_bowtie.log"

    bowtie \
      -f -v 0 -p 2 -a --best --strata \
      ref/intron \
      "reads_mapping/${sample}.mapped.fa" \
      > "reads_mapping/${sample}.intron.bwt" \
      2> "logs/${sample}.intron_bowtie.log"

    bowtie \
      -f -v 0 -p 2 -a --best --strata \
      ref/exon \
      "reads_mapping/${sample}.mapped.fa" \
      > "reads_mapping/${sample}.exon.bwt" \
      2> "logs/${sample}.exon_bowtie.log"
done
