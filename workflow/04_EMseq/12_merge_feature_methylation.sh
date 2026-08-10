#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

HELPER="software/batmeth2merge.R"
[[ -f "$HELPER" ]] || { echo "[ERROR] Missing $HELPER" >&2; exit 1; }

pushd FeatureML >/dev/null

: > genes.AverMethylevel.txt
: > TEs.AverMethylevel.txt
: > genes.Promoter.c.txt

for sample in "${SAMPLES[@]}"; do
    echo "${sample}.genes.AverMethylevel.txt" >> genes.AverMethylevel.txt
    echo "${sample}.TEs.AverMethylevel.txt" >> TEs.AverMethylevel.txt
    echo "${sample}.genes.Promoter.c.txt" >> genes.Promoter.c.txt
done

Rscript ../"$HELPER" -p genes.AverMethylevel.txt -s 1 -e 6
Rscript ../"$HELPER" -p TEs.AverMethylevel.txt -s 1 -e 6
Rscript ../"$HELPER" -p genes.Promoter.c.txt -s 1 -e 6

popd >/dev/null
