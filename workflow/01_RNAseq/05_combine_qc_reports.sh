#!/usr/bin/env bash
set -euo pipefail

# These are the exact helper-script calls recorded in the supplied analysis log.
# Copy the original helper scripts into ./scripts/ before running.

Rscript scripts/combine_fastp_report.R \
    -i clean_data/ \
    -o clean_data/

python scripts/combine_hisat2_report.py \
    -i aligned/ \
    -o aligned/
