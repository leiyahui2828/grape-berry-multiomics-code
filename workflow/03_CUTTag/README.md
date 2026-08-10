# CUT&Tag histone-modification analysis code

This directory is a cleaned, non-Snakemake organization of the CUT&Tag
commands represented by the supplied `cut_tag_code` dry-run and generic
`Snakefile`.

The attached dry-run is specifically for **H3K27ac** (five developmental
stages × two biological replicates). The same repository structure is prepared
for the other three histone marks used in the study:

| Histone mark | Peak class |
|---|---|
| H3K27ac | narrowPeak |
| H3K4me3 | narrowPeak |
| H3K9ac | narrowPeak |
| H3K27me3 | broadPeak |

The narrow/broad distinction is handled explicitly by the peak-calling,
peak-filtering, reproducibility and sequence-extraction scripts.

## Directory structure

```text
03_CUTTag/
├── README.md
├── METHOD_AUDIT.md
├── 00_samples.tsv
├── 00_config.sh
├── 01_fastp.sh
├── 02_bowtie2_build.sh
├── 03_bowtie2_align.sh
├── 04_samtools_sort.sh
├── 05_filter_bam.sh
├── 06_qc.sh
├── 07_bigwig.sh
├── 08_deeptools_profiles.sh
├── 09_MACS3_peak_calling.sh
├── 10_filter_peaks_qscore.sh
├── 11_select_reproducible_peaks.sh
├── 12_peak_annotation_uropa.sh
├── 13_merge_consensus_peaks.sh
├── 14_featureCounts.sh
├── 15_build_count_matrix.sh
├── 16_edgeR_differential_peaks.sh
├── 17_extract_peak_sequence.sh
├── run_mark.sh
├── scripts/
│   └── README.md
├── software/
│   └── RunFeatureCounts/
│       └── README.md
└── ref/
    └── README.md
```

## Scope of the organized workflow

```text
FASTQ
  ↓
fastp
  ↓
Bowtie2
  ↓
SAMtools
  ↓
alignmentSieve
  ↓
filtered BAM
  ├── bamCoverage / deepTools QC
  ↓
MACS3 peak calling
  ├── narrowPeak: H3K27ac / H3K4me3 / H3K9ac
  └── broadPeak:  H3K27me3
  ↓
peak q-score filtering
  ↓
biological-replicate reproducibility
  ↓
stage-level peaks
  ↓
UROPA annotation
  ↓
merged consensus peak set
  ↓
featureCounts
  ↓
peak count matrix
  ↓
edgeR differential peak analysis
```

## H3K27ac parameters directly supported by the uploaded dry-run

```text
fastp:
  --thread 3
  --qualified_quality_phred 20
  --length_required 50

Bowtie2:
  -p 4
  --no-mixed
  --no-discordant

alignmentSieve:
  --ignoreDuplicates
  --minMappingQuality 25
  --samFlagExclude 260

bamCoverage:
  --effectiveGenomeSize 494873210
  --binSize 10
  --normalizeUsing RPKM

MACS3:
  --format BAM
  --keep-dup all
  --qvalue 0.075
  --cutoff-analysis
  --nomodel
  --shift -100
  --extsize 200
  --SPMR
  -B

post-MACS3:
  q-score column 9 > 3

H3K27ac replicate selection:
  IDR threshold = 0.05
  rank = p.value

Differential peaks:
  edgeR via Trinity run_DE_analysis.pl
  dispersion = 0.1
```

## Run one histone mark

For the attached H3K27ac workflow:

```bash
bash run_mark.sh H3K27ac
```

All outputs are separated by mark:

```text
results/H3K27ac/
results/H3K4me3/
results/H3K9ac/
results/H3K27me3/
```

This prevents peaks from different histone modifications from being merged
into the same consensus peak set.

## Other three histone marks

The directory already knows the expected peak class:

```text
H3K4me3  -> narrowPeak
H3K9ac   -> narrowPeak
H3K27me3 -> broadPeak
```

However, the uploaded dry-run does **not** contain their exact historical
MACS3-extra parameters or peak-selection setting. Before running those marks,
restore/check their original config files and provide the actual MACS3 values
to `09_MACS3_peak_calling.sh`.

See `METHOD_AUDIT.md` for the items that must be reconciled before the
GitHub/Zenodo `v1.0.0` release.
