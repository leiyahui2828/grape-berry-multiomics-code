#!/usr/bin/env bash

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.tsv}"

# Reference/database locations used by the organized workflow.
GENOME="${GENOME:-ref/genome.fa}"
GENOME_GFF="${GENOME_GFF:-ref/genome.gff3}"
REPEAT_GFF="${REPEAT_GFF:-ref/repeat.gff3}"

MIRBASE_DIR="${MIRBASE_DIR:-ref/database/miRbase}"
RFAM_DIR="${RFAM_DIR:-ref/database/rfam}"

SPECIES_PREFIX="${SPECIES_PREFIX:-vvi}"

# Historical local software entry points.
MAPPER_PL="${MAPPER_PL:-software/mirdeep2-0.1.3/bin/mapper.pl}"
QUANTIFIER_PL="${QUANTIFIER_PL:-software/mirdeep2-0.1.3/bin/quantifier.pl}"
MIRDP2_PIPELINE="${MIRDP2_PIPELINE:-software/miRDP2-v1.1.4/miRDP2-v1.1.4_pipeline.bash}"
PSROBOT="${PSROBOT:-software/psRobot_v1.2/bin/psRobot_tar}"
TARGETFINDER="${TARGETFINDER:-software/targetfinder_threads.pl}"

get_samples() {
    awk -F '\t' 'NR>1 {print $1}' "$SAMPLES_FILE"
}

get_stage() {
    awk -F '\t' -v s="$1" 'NR>1 && $1==s {print $2}' "$SAMPLES_FILE"
}
