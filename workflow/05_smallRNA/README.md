# small RNA / miRNA analysis

This directory organizes the commands in the supplied
`5.small_rna组学数据分析(2).txt` as ordinary Bash scripts.

Study design:

```text
PN4, PN7, PN10, PN15, PN18
2 biological replicates per stage
10 single-end small-RNA libraries
```

## Workflow overview

```text
miRBase / Rfam / PN40024 reference preparation
        ↓
FastQC
        ↓
Cutadapt
        ↓
Trim Galore (18–30 nt, Q30)
        ↓
miRDeep2 mapper.pl
        ↓
collapsed non-redundant FASTA
        ↓
Bowtie classification
  ├── genome
  ├── known grape miRNA precursor
  ├── Rfam non-miRNA ncRNA
  ├── repeat
  ├── intron
  └── exon/CDS-like reference
        ↓
srna_anno.pl
        ↓
known miRNA + candidate novel-miRNA reads
        ↓
miRDP2 novel miRNA prediction
        ↓
miRNA sequence/base-content analyses
        ↓
miRNA family / phylogeny / motif analyses
        ↓
miRDeep2 quantifier.pl
        ↓
miRNA count / TPM matrices
        ↓
DESeq2 differential miRNA analysis
        ↓
Target prediction
  ├── psRobot
  ├── TargetFinder
  └── psRNATarget
        ↓
intersection of supported target genes
```

## Directory structure

```text
05_smallRNA/
├── README.md
├── METHOD_AUDIT.md
├── 00_samples.tsv
├── 00_config.sh
├── 01_prepare_miRBase.sh
├── 02_prepare_Rfam.sh
├── 03_prepare_reference.sh
├── 04_fastqc_trim.sh
├── 05_collapse_reads.sh
├── 06_length_statistics.sh
├── 07_map_genome.sh
├── 08_map_annotation_databases.sh
├── 09_merge_sRNA_annotation.sh
├── 10_prepare_novel_miRNA_reads.sh
├── 11_miRDP2_predict.sh
├── 12_rename_novel_miRNA.sh
├── 13_miRNA_base_content.sh
├── 14_miRNA_family_annotation.sh
├── 15_phylogeny_and_motif.sh
├── 16_miRNA_quantification.sh
├── 17_expression_plots.sh
├── 18_DESeq2_miRNA.sh
├── 19_target_prediction.sh
├── 20_target_intersection.sh
├── run_core_upstream.sh
├── scripts/
│   └── README.md
├── software/
│   └── README.md
└── ref/
    └── README.md
```

## Key historical parameters

### Read filtering

```text
Cutadapt:
  -g ACGATC
  -a AGATCG

Trim Galore:
  --small_rna
  --length 18
  --max_length 30
  --quality 30
  --stringency 3
  --phred33
  --cores 6
  --dont_gzip
```

The source file first ran Trim Galore without `--max_length 30`, explicitly
noted that this failed to remove reads >30 nt, and then re-ran the command with
`--max_length 30`. The organized workflow uses the corrected final command.

### Genome and annotation mapping

```text
Genome:
  Bowtie -f -v 0 -k 1

Known grape miRNA precursor:
  Bowtie -f -v 0 -p 8 -a --best --strata

Rfam non-miRNA ncRNA:
  Bowtie -f -v 1 -p 8 -a --best --strata

Repeat / intron / exon:
  Bowtie -f -v 0 -p 2 -a --best --strata
```

### Novel miRNA prediction

```text
miRDP2 pipeline v1.1.4
--genome ref/genome.fa
--index ref/genome
--fasta
```

### miRNA quantification

```text
miRDeep2 quantifier.pl
-g 0
```

### Differential miRNA expression

```text
DESeq2 via script/run_DE_analysis.pl
```

### Target prediction

```text
psRobot
TargetFinder
psRNATarget (online result imported locally)
```

## Important

Several custom helper scripts referenced by the original analysis were not
included in the uploaded text file. Restore the exact historical versions listed
in `scripts/README.md` before the final GitHub/Zenodo release.

Also read `METHOD_AUDIT.md`: the source code contains several copy/paste and
Bowtie/Bowtie2-index inconsistencies that should be verified before tagging the
repository as `v1.0.0`.
