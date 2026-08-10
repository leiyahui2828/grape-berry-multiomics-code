#!/usr/bin/env bash
set -euo pipefail

mkdir -p logs

hicFindRestSite \
    --fasta ref/genome.fa \
    --searchPattern GATC \
    -o ref/rest_site.bed \
    > logs/hicFindRestSite.log 2>&1
