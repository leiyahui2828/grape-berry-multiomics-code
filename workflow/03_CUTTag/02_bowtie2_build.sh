#!/usr/bin/env bash
set -euo pipefail

source ./00_config.sh
mkdir -p ref logs

bowtie2-build \
    --threads 6 \
    "$GENOME" \
    ref/genome \
    > logs/bowtie2_build.log 2>&1
