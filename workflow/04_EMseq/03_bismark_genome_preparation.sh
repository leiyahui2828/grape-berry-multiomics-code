#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mkdir -p logs
run_meth_tools bismark_genome_preparation ref/ \
  > logs/bismark_genome_preparation.log 2>&1
