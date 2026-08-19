#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mkdir -p ref logs

[[ -f "$GENOME" ]] || { echo "[ERROR] Missing $GENOME" >&2; exit 1; }
[[ -f "$GTF" ]] || { echo "[ERROR] Missing $GTF" >&2; exit 1; }

gffread "$GTF" -g "$GENOME" -w "$TRANSCRIPTS"

gffread -C --keep-genes "$GTF" -T > ref/gene.coding.gff3

agat_sp_keep_longest_isoform.pl \
  --cpu 10 \
  --gff ref/gene.coding.gff3 \
  -o "$LONGEST_GFF"

agat_sp_extract_sequences.pl \
  --gff "$LONGEST_GFF" \
  --fasta "$GENOME" \
  --protein \
  -o ref/genome.protein.fa

awk -F '\t' '$3=="transcript"{
  match($9,/ID=([^;]+)/,a);
  match($9,/Parent=([^;]+)/,b);
  print a[1]"\t"b[1]
}' "$LONGEST_GFF" > ref/gene.final.mapid

seqkit replace \
  -k ref/gene.final.mapid \
  -p '^(\S+)' \
  -r '{kv}' \
  ref/genome.protein.fa \
  > ref/genome.protein_geneID.fa

agat_sp_extract_sequences.pl \
  --gff "$LONGEST_GFF" \
  --fasta "$GENOME" \
  --type five_prime_utr \
  -o ref/genome.five_prime_utr.fa

seqkit replace \
  -k ref/gene.final.mapid \
  -p '^(\S+)' \
  -r '{kv}' \
  ref/genome.five_prime_utr.fa \
  > ref/genome.five_prime_utr_geneID.fa

if [[ -f scripts/fasta_GC_NMFE.py ]]; then
  python scripts/fasta_GC_NMFE.py \
    --fasta ref/genome.five_prime_utr_geneID.fa \
    --out ref/genome.five_prime_utr_geneID.GC_NMFE.tab \
    --threads 20
fi
