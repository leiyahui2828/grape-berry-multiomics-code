#!/usr/bin/env bash
set -euo pipefail

head -1 DM_Analysis/CpG_region/PN7_vs_PN4.dmr \
  > DM_Analysis/merged.dmr

cat DM_Analysis/C*_region/*.dmr \
  | grep -v "start" \
  >> DM_Analysis/merged.dmr
