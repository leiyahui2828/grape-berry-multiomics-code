#!/usr/bin/env bash
set -euo pipefail

mkdir -p logs

bwa index \
    -p ref/genome \
    ref/genome.fa \
    > logs/bwa_index.log 2>&1
