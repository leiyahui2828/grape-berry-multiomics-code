#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh

RNA_TPM="${RNA_TPM:-../01_RNAseq/quanti/genes.TMM.TPM.matrix}"
RNA_COUNTS="${RNA_COUNTS:-../01_RNAseq/quanti/genes.counts.matrix}"
RNA_DE="${RNA_DE:-../01_RNAseq/DE_analysis/de_result.tsv}"

DIR="$OUTDIR/RNA_Ribo_integration"
mkdir -p "$DIR"

for f in "$RNA_TPM" "$RNA_COUNTS" \
         "$OUTDIR/RPF_quant/matrix.TPM.not_cross_norm" \
         "$OUTDIR/RPF_quant/matrix.counts.matrix"; do
    [[ -f "$f" ]] || { echo "[ERROR] Missing $f" >&2; exit 1; }
done

sed '1s/\t/\tRIBO_/g' \
  "$OUTDIR/RPF_quant/matrix.TPM.not_cross_norm" \
  > "$DIR/Ribo-seq.TPM.matrix"

sed '1s/\t/\tRNA_/g' \
  "$RNA_TPM" \
  > "$DIR/RNA-seq.TPM.matrix"

sed '1s/^/GeneID/' "$OUTDIR/RPF_quant/matrix.counts.matrix" \
  | sed '1s/\t/\tRIBO_/g' \
  > "$DIR/Ribo-seq.counts.matrix"

sed '1s/^/GeneID/' "$RNA_COUNTS" \
  | sed '1s/\t/\tRNA_/g' \
  > "$DIR/RNA-seq.counts.matrix"

echo -e "SampleID\tCondition\tSeqType" > "$DIR/sample_info_all.txt"

for seqtype in RIBO RNA; do
  for stage in PN4 PN7 PN10 PN15 PN18; do
    for rep in 1 2 3; do
      printf "%s_%s_rep%s\t%s\t%s\n" \
        "$seqtype" "$stage" "$rep" "$stage" "$seqtype" \
        >> "$DIR/sample_info_all.txt"
    done
  done
done

[[ -f "$RNA_DE" ]] && ln -sf "$(realpath "$RNA_DE")" "$DIR/RNA-seq.DE_results"
ln -sf "$(realpath "$OUTDIR/RPF_quant/all.RPF.DE_results")" "$DIR/Ribo-seq.DE_results"
