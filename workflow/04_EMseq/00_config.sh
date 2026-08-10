#!/usr/bin/env bash
SAMPLES_FILE="${SAMPLES_FILE:-00_samples.tsv}"
METH_TOOLS_SIF="${METH_TOOLS_SIF:-sifs/methtools20220515.sif}"
GENOME="${GENOME:-ref/genome.fa}"
GFF="${GFF:-ref/genes.gff}"
GTF="${GTF:-ref/genes.gtf}"
TE_GFF="${TE_GFF:-ref/TE.gff}"

FASTP_THREADS="${FASTP_THREADS:-5}"
BISMARK_THREADS="${BISMARK_THREADS:-4}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-4}"

METH_BASEQ="${METH_BASEQ:-20}"
METH_COVERAGE="${METH_COVERAGE:-4}"
METH_NC="${METH_NC:-1}"
METH_REGION="${METH_REGION:-600}"
METH_STEP="${METH_STEP:-50000}"

FEATURE_DISTANCE="${FEATURE_DISTANCE:-2000}"
FEATURE_STEP="${FEATURE_STEP:-0.01}"

run_meth_tools() {
    singularity exec "$METH_TOOLS_SIF" "$@"
}

get_samples() {
    awk -F '\t' 'NR>1 {print $1}' "$SAMPLES_FILE"
}
