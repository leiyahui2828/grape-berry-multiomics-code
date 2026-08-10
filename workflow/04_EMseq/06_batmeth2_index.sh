#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mkdir -p logs
run_meth_tools BatMeth2 index -g "$GENOME" \
  > logs/genome_BatMeth_index.log 2>&1
