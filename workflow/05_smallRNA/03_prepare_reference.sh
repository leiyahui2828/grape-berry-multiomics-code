#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh

mkdir -p ref

for f in "$GENOME" "$GENOME_GFF" "$REPEAT_GFF"; do
    [[ -f "$f" ]] || { echo "[ERROR] Missing reference file: $f" >&2; exit 1; }
done

# Genome index used in historical Bowtie mapping and miRDP2.
bowtie-build "$GENOME" ref/genome

# Repeat reference.
awk '{print $1 "\t" $4 "\t" $5}' "$REPEAT_GFF" > ref/repeat.bed
bedtools getfasta -fi "$GENOME" -fo ref/repeat.fa.tmp -bed ref/repeat.bed

seqtk comp ref/repeat.fa.tmp \
  | awk '$3+$4+$5+$6>0 && $2>20 {print $1}' \
  > ref/repeat.fa.tmp.id

seqtk subseq ref/repeat.fa.tmp ref/repeat.fa.tmp.id > ref/repeat.fa

# Historical source used bowtie2-build here.
bowtie2-build ref/repeat.fa ref/repeat

# Intron reference.
[[ -f scripts/get_intron.pl ]] || {
    echo "[ERROR] Missing scripts/get_intron.pl" >&2
    exit 1
}

perl scripts/get_intron.pl "$GENOME_GFF" > ref/intron.bed
bedtools getfasta -fi "$GENOME" -fo ref/intron.fa -bed ref/intron.bed -name -s

# Historical source used bowtie2-build here.
bowtie2-build ref/intron.fa ref/intron

# Exon/CDS-like reference:
# The supplied source explicitly states that this step was uncertain.
# Do not invent its construction. If ref/exon.fa exists, preserve the
# historical Bowtie2-index command.
if [[ -f ref/exon.fa ]]; then
    bowtie2-build ref/exon.fa ref/exon
else
    echo "[WARNING] ref/exon.fa not found; exon index not built." >&2
    echo "[WARNING] See METHOD_AUDIT.md before publication release." >&2
fi
