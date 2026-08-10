#!/usr/bin/env bash
set -euo pipefail

mkdir -p logs

hisat2-build \
    ref/genome.fa \
    genome.hisat2 \
    > logs/build_hisat2_index.log 2>&1
