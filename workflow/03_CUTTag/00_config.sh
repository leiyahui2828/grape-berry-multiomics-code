#!/usr/bin/env bash

# Common parameters reconstructed from the supplied H3K27ac dry-run.
#
# IMPORTANT:
# The uploaded execution log contains H3K27ac only. Common upstream parameters
# below are therefore source-supported for H3K27ac. Exact mark-specific MACS3
# extras / reproducibility settings for H3K4me3, H3K9ac and H3K27me3 should be
# checked against their historical config files before publication release.

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.tsv}"

GENOME="${GENOME:-ref/genome.fa}"
GTF="${GTF:-ref/genes_longest.gtf}"
TSS_TES_BED="${TSS_TES_BED:-ref/genes_longest.bed6}"
TSS_TES_SHUFFLE_BED="${TSS_TES_SHUFFLE_BED:-ref/genes_shuffle.bed6}"

EFFECTIVE_GENOME_SIZE="${EFFECTIVE_GENOME_SIZE:-494873210}"

FASTP_THREADS="${FASTP_THREADS:-3}"
BOWTIE2_THREADS="${BOWTIE2_THREADS:-4}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-4}"
FILTER_THREADS="${FILTER_THREADS:-4}"
BAMCOVERAGE_THREADS="${BAMCOVERAGE_THREADS:-2}"

MIN_MAPQ="${MIN_MAPQ:-25}"
QSCORE_CUTOFF="${QSCORE_CUTOFF:-3}"
IDR_THRESHOLD="${IDR_THRESHOLD:-0.05}"

# The historical H3K27ac dry-run used:
#   --blackListFileName ref/hg38.blacklist.bed
#
# Because the study organism is grape, this filename must be verified before
# reuse. The publication-oriented scripts therefore do NOT enable it by default.
#
# If the file is confirmed to be a grape blacklist that was merely misnamed:
#   export BLACKLIST=ref/hg38.blacklist.bed
#
BLACKLIST="${BLACKLIST:-}"

# H3K27ac source-supported MACS3 settings from the attached dry-run.
H3K27AC_MACS3_FORMAT="${H3K27AC_MACS3_FORMAT:-BAM}"
H3K27AC_MACS3_EXTRA="${H3K27AC_MACS3_EXTRA:---nomodel --shift -100 --extsize 200 --SPMR -B}"

get_peak_type() {
    case "$1" in
        H3K27ac|H3K4me3|H3K9ac) echo "narrow" ;;
        H3K27me3) echo "broad" ;;
        *)
            echo "[ERROR] Unsupported histone mark: $1" >&2
            return 1
            ;;
    esac
}

get_samples() {
    local mark="$1"
    awk -F '\t' -v m="$mark" 'NR>1 && $2==m {print $1}' "$SAMPLES_FILE"
}

get_stages() {
    printf '%s\n' PN4 PN7 PN10 PN15 PN18
}

result_dir() {
    echo "results/$1"
}
