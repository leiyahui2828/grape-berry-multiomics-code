#!/usr/bin/env bash

# Grape-adapted Ribo-seq workflow.
# Defaults below are inherited from the uploaded Arabidopsis reference workflow
# and must be verified against the actual grape library preparation / final code.

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.tsv}"
SAMPLE_INFO="${SAMPLE_INFO:-sample_info.txt}"
CONTRASTS_FILE="${CONTRASTS_FILE:-contrasts.txt}"

RAW_DIR="${RAW_DIR:-raw_data}"
OUTDIR="${OUTDIR:-results}"

GENOME="${GENOME:-ref/genome.fa}"
GTF="${GTF:-ref/genome.gtf}"
TRANSCRIPTS="${TRANSCRIPTS:-ref/transcripts.fa}"
LONGEST_GFF="${LONGEST_GFF:-ref/genome.longest.gff3}"

NCRNA_FASTA="${NCRNA_FASTA:-ref/ncRNA.fasta}"
NCRNA_INDEX="${NCRNA_INDEX:-ref/ncRNA}"
STAR_INDEX="${STAR_INDEX:-ref/STAR_index}"

RIBO_ADAPTER="${RIBO_ADAPTER:-TGGAATTCTCGGGTGCCAA}"
FASTP_MIN_LEN="${FASTP_MIN_LEN:-20}"
FASTP_MAX_LEN="${FASTP_MAX_LEN:-40}"
FASTP_THREADS="${FASTP_THREADS:-4}"

NCRNA_MISMATCH="${NCRNA_MISMATCH:-2}"
BOWTIE_THREADS="${BOWTIE_THREADS:-4}"

STAR_THREADS="${STAR_THREADS:-8}"
STAR_MULTIMAP_NMAX="${STAR_MULTIMAP_NMAX:-1}"
STAR_MISMATCH_NMAX="${STAR_MISMATCH_NMAX:-2}"

RPF_MIN_LEN="${RPF_MIN_LEN:-26}"
RPF_MAX_LEN="${RPF_MAX_LEN:-38}"

RPF_STRAND="${RPF_STRAND:-1}"
RPF_FEATURE="${RPF_FEATURE:-CDS}"
RPF_ATTRIBUTE="${RPF_ATTRIBUTE:-gene_id}"
RPF_IS_PAIRED="${RPF_IS_PAIRED:-FALSE}"

get_samples() {
    awk -F '\t' 'NR>1 {print $1}' "$SAMPLES_FILE"
}
get_fastq() {
    awk -F '\t' -v s="$1" 'NR>1 && $1==s {print $4}' "$SAMPLES_FILE"
}
