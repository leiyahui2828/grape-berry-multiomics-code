# 一.安装

1. 安装依赖R包 argparser, Rsubread, limma, edgeR
2. 安装本软件

# 二.使用

对 GTF 中的 GENE 定量
```
Rscript RunFeatureCounts/run-featurecounts.R -b ../21.Mapping/BLO_S1_LD1.bam -g ../12.ref/genes.gtf -f exon -a gene_id -i FALSE -s 2 -o BLO_S1_LD1
```
对 GTF 中的 ISOFORM 定量
```
Rscript RunFeatureCounts/run-featurecounts.R -b ../21.Mapping/BLO_S1_LD1.bam -g ../12.ref/genes.gtf -f exon -a transcript_id -i FALSE -s 2 -o BLO_S1_LD1
```


对 SAF 定量
```
Rscript RunFeatureCounts/run-featurecounts.R -b aligned/FoxA1_E2_rep1_sorted_filtered.bam --saf merged_peaks.saf --isPairedEnd FALSE -o FoxA1_E2_rep1
```
