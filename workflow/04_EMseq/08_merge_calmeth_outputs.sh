#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

HELPER="software/batmeth2merge.R"
[[ -f "$HELPER" ]] || { echo "[ERROR] Missing $HELPER" >&2; exit 1; }

mkdir -p methylation_merged logs

make_list () {
    suffix="$1"
    outfile="$2"
    : > "$outfile"
    for sample in "${SAMPLES[@]}"; do
        f="site_samll_window_ML/${sample}.${suffix}"
        [[ -f "$f" ]] || { echo "[ERROR] Missing $f" >&2; exit 1; }
        echo "$f" >> "$outfile"
    done
}

make_list "methBins.txt" methylation_merged/methBins.txt
make_list "Region.CG.txt" methylation_merged/Region.CG.txt
make_list "Region.CHG.txt" methylation_merged/Region.CHG.txt
make_list "Region.CHH.txt" methylation_merged/Region.CHH.txt

pushd methylation_merged >/dev/null

Rscript ../"$HELPER" -p methBins.txt -s 1 -e 9
sed -i 's/\.//g' merged.methBins.txt

for context in CG CHG CHH; do
    Rscript ../"$HELPER" -p "Region.${context}.txt" -s 1 -e 9 \
      > "../logs/batmeth2merge_Region.${context}.log" 2>&1
    sed -i 's/_$//' "merged.Region.${context}.txt"
done

cat merged.Region.*.txt > merged.Region.C.txt
popd >/dev/null
