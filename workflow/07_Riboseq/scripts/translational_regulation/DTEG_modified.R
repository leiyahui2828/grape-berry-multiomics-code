suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(edgeR)
  library(argparser, quietly = TRUE)
})

## Input files
# Calculating differential translation genes requires the count matrices from Ribo-seq and RNA-seq.
# These should be raw counts, not normalized values.
# Sample information file must contain sample, Condition, SeqType.
# If batch == 1, it must also contain Batch.


### Get command line arguments using argparser

p <- arg_parser(
  description = paste(
    "Differential translation efficiency analysis using DESeq2.",
    "deltaTE log2FoldChange = (RIBO treat vs control) - (RNA treat vs control)"
  )
)

p <- add_argument(
  p,
  "--ribo_count",
  help = "Ribo-seq raw count matrix. The first column must be gene ID.",
  type = "character",
  default = NA
)

p <- add_argument(
  p,
  "--rna_count",
  help = "RNA-seq raw count matrix. The first column must be gene ID.",
  type = "character",
  default = NA
)

p <- add_argument(
  p,
  "--sample_info",
  help = "Sample information file. Required columns: sample, Condition, SeqType. If --batch 1, Batch is also required.",
  type = "character",
  default = NA
)

p <- add_argument(
  p,
  "--outdir",
  help = "Output directory.",
  type = "character",
  default = "Results"
)

p <- add_argument(
  p,
  "--batch",
  help = "Whether to include Batch in the model: 0 or 1.",
  type = "integer",
  default = 0
)

p <- add_argument(
  p,
  "--control",
  help = "Control condition name.",
  type = "character",
  default = NA
)

p <- add_argument(
  p,
  "--treat",
  help = "Treatment condition name.",
  type = "character",
  default = NA
)

p <- add_argument(
  p,
  "--padj_cutoff",
  help = "Adjusted P-value cutoff.",
  type = "numeric",
  default = 0.05
)

p <- add_argument(
  p,
  "--log2fc_cutoff",
  help = "Absolute log2 fold change cutoff. Set 0 to disable log2FC filtering.",
  type = "numeric",
  default = 0
)

argv <- parse_args(p)

# If Rscript deltaTE_optimized.R is run without arguments,
# argparser will print help automatically only when required arguments are defined without default.
# Here we use default = NA for compatibility, so manually show help when key arguments are missing.

if (
  is.na(argv$ribo_count) ||
  is.na(argv$rna_count) ||
  is.na(argv$sample_info) ||
  is.na(argv$control) ||
  is.na(argv$treat)
) {
  print(p)
  quit(save = "no", status = 1)
}

# Input filenames and parameters

ribo_file <- argv$ribo_count
rna_file <- argv$rna_count
data_file <- argv$sample_info
outdir <- argv$outdir
batch <- argv$batch
control <- argv$control
treat <- argv$treat
padj_cutoff <- argv$padj_cutoff
log2fc_cutoff <- argv$log2fc_cutoff

if (!batch %in% c(0, 1)) {
  stop("Batch presence should be indicated by 0 or 1 only.", call. = FALSE)
}

if (is.na(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff >= 1) {
  stop("padj_cutoff should be a number between 0 and 1.", call. = FALSE)
}

if (is.na(log2fc_cutoff) || log2fc_cutoff < 0) {
  stop("log2fc_cutoff should be a non-negative number.", call. = FALSE)
}

if (control == treat) {
  stop("control and treat should be different.", call. = FALSE)
}

## 增加差异倍数参数
sig_filter <- function(res_obj) {
  sig = !is.na(res_obj$padj) & res_obj$padj < padj_cutoff
  
  if (log2fc_cutoff > 0) {
    sig = sig & !is.na(res_obj$log2FoldChange) &
      abs(res_obj$log2FoldChange) >= log2fc_cutoff
  }
  
  sig
}

## 增加输出函数
# Add GeneID, control, treat and comparison columns before writing DESeq2 results
write_deseq_result <- function(res_obj, outfile, control, treat) {
  
  res_df <- as.data.frame(res_obj)
  
  res_df <- data.frame(
    GeneID = rownames(res_df),
    control = control,
    treat = treat,
    comparison = paste0(treat, "_vs_", control),
    res_df,
    check.names = FALSE
  )
  
  write.table(
    res_df,
    outfile,
    quote = FALSE,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
  )
}


### Read count matrices
# The first column is gene ID.
# check.names = FALSE is used to preserve sample names containing "-".
### Read count matrices and sample information

# 读取 count 文件：第一列为 gene_id，后续列为样品 raw count
read_count <- function(file, name = "count") {
  x = read.delim(
    file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (ncol(x) < 1) {
    stop(paste0(name, " file should contain gene ID column and sample count columns."),call. = FALSE)
  }
  
  if (any(duplicated(rownames(x)))) {
    stop(paste0("Duplicated gene IDs found in ", name, " file."), call. = FALSE )
  }
  
  if (any(is.na(rownames(x))) || any(rownames(x) == "")) {
    stop(paste0(name, " file contains empty gene IDs."),call. = FALSE )
  }
  
  x = as.matrix(x)
  suppressWarnings(storage.mode(x) <- "numeric")
  
  if (any(is.na(x))) {
    stop( paste0(name, " matrix contains NA or non-numeric values."), call. = FALSE )
  }
  
  if (any(x < 0)) {
    stop( paste0(name, " matrix contains negative values."),call. = FALSE )
  }
  
  if (any(abs(x - round(x)) > 1e-6)) {
    warning(paste0( name," matrix contains non-integer values, rounded to integers for DESeq2 analysis."),call. = FALSE
    )
  }
  
  round(x)
}

ribo = read_count(ribo_file, "Ribo-seq count")
rna  = read_count(rna_file, "RNA-seq count")

# 过滤地表达基因
# ribo = ribo[rowSums(cpm(ribo) > 1) >= 2,]
# rna = rna[rowSums(cpm(rna) > 1) >= 2,]


# 按 gene_id 对齐 Ribo-seq 和 RNA-seq
common_genes = intersect(rownames(ribo), rownames(rna))
ribo_only_genes = setdiff(rownames(ribo), rownames(rna))
rna_only_genes = setdiff(rownames(rna), rownames(ribo))

if (length(common_genes) == 0) {
  stop("No common genes found between Ribo-seq and RNA-seq count files.", call. = FALSE)
}

ribo = ribo[common_genes, , drop = FALSE]
rna  = rna[common_genes, , drop = FALSE]

merge = cbind(ribo, rna)

if (any(duplicated(colnames(merge)))) {
  stop("Duplicated sample names found after merging Ribo-seq and RNA-seq count matrices.", call. = FALSE)
}

# 读取样品信息表
coldata = read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)

required_cols = c("SampleID", "Condition", "SeqType")
if (batch == 1) {
  required_cols = c(required_cols, "Batch")
}

missing_cols = setdiff(required_cols, colnames(coldata))
if (length(missing_cols) > 0) {
  stop(
    paste0("Sample information file is missing required columns: ",paste(missing_cols, collapse = ", ")),
    call. = FALSE
  )
}

if (any(duplicated(coldata$SampleID))) {
  stop("Duplicated sample names found in sample information file.", call. = FALSE)
}

if (!all(c(control, treat) %in% coldata$Condition)) {
  stop("control or treat condition was not found in sample information file.", call. = FALSE)
}

if (!all(coldata$SeqType %in% c("RNA", "RIBO"))) {
  stop("SeqType must only contain RNA and RIBO.", call. = FALSE)
}

if (batch == 1 && any(is.na(coldata$Batch) | coldata$Batch == "")) {
  stop("Batch column contains NA or empty values.", call. = FALSE)
}

rownames(coldata) = coldata$SampleID

# sample_info 可以包含多组 Condition，这里只保留本次比较的 control 和 treat
coldata = coldata[coldata$Condition %in% c(control, treat), , drop = FALSE]

# 只保留 count 矩阵和 sample_info 中都存在的样品
keep_samples = intersect(colnames(merge), rownames(coldata))

if (length(keep_samples) == 0) {
  stop("No matched samples found between count matrix and selected conditions.", call. = FALSE)
}

merge = merge[, keep_samples, drop = FALSE]
coldata = coldata[keep_samples, , drop = FALSE]

if (!identical(colnames(merge), rownames(coldata))) {
  stop("Sample names in count matrix and sample information file are not matched.", call. = FALSE)
}

# 同步更新单独的 ribo 和 rna 矩阵，供后续单独差异分析和 PCA 使用
ribo = ribo[, intersect(colnames(ribo), keep_samples), drop = FALSE]
rna  = rna[, intersect(colnames(rna), keep_samples), drop = FALSE]

# 检查每个 SeqType × Condition 是否都有样品
sample_table = table(coldata$SeqType, coldata$Condition)
print(sample_table)

if (!all(c("RNA", "RIBO") %in% rownames(sample_table))) {
  stop("Both RNA and RIBO samples are required after condition filtering.", call. = FALSE)
}

if (!all(c(control, treat) %in% colnames(sample_table))) {
  stop("Both control and treat conditions are required after condition filtering.", call. = FALSE)
}

if (any(sample_table[c("RNA", "RIBO"), c(control, treat)] == 0)) {
  stop("Each SeqType and Condition combination must contain at least one sample.", call. = FALSE)
}

if (any(sample_table[c("RNA", "RIBO"), c(control, treat)] < 2)) {
  warning("Some SeqType and Condition combinations have fewer than 2 samples. DESeq2 results may be unreliable.")
}

# 保存原始 Condition 名称用于绘图
coldata$Condition_original = factor(
  as.character(coldata$Condition),
  levels = c(as.character(control), as.character(treat))
)

# 将 Condition 标准化为 control/treat，避免条件名中含有特殊字符导致系数名不稳定
coldata$Condition = ifelse(coldata$Condition == control, "control", "treat")
coldata$Condition = factor(coldata$Condition, levels = c("control", "treat"))

# SeqType 固定为 RNA/RIBO，RNA 作为参考水平
coldata$SeqType = factor(coldata$SeqType, levels = c("RNA", "RIBO"))
coldata$SeqType = relevel(coldata$SeqType, "RNA")

if (batch == 1) {
  coldata$Batch = factor(coldata$Batch)
}

#head(merge)
head(coldata)

## Detecting differential translation regulation
### DESeq2 object with batch and interaction term in the design

if (batch == 1) {
  ddsMat = DESeqDataSetFromMatrix(
    countData = merge,
    colData = coldata,
    design = ~ Batch + Condition + SeqType + Condition:SeqType
  )
} else if (batch == 0) {
  ddsMat = DESeqDataSetFromMatrix(
    countData = merge,
    colData = coldata,
    design = ~ Condition + SeqType + Condition:SeqType
  )
} else {
  stop("Batch presence should be indicated by 0 or 1 only", call. = FALSE)
}

ddsMat = DESeq(ddsMat)
res_names = resultsNames(ddsMat)
print(res_names)

# After Condition is standardized as control/treat and SeqType is RNA/RIBO,
# the interaction coefficient should usually be Conditiontreat.SeqTypeRIBO.
interaction_term = grep("^Conditiontreat\\.SeqTypeRIBO$", res_names, value = TRUE)

if (length(interaction_term) != 1) {
  interaction_term = grep("Condition.*treat.*SeqType.*RIBO", res_names, value = TRUE)
}

if (length(interaction_term) != 1) {
  stop(
    paste0(
      "Cannot identify the Condition:SeqType interaction term.\n",
      "Available resultsNames are:\n",
      paste(res_names, collapse = "\n")
    ),
    call. = FALSE
  )
}

# Output directories
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "fold_changes"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "gene_lists"), showWarnings = FALSE, recursive = TRUE)

# Input check information
check_file = file.path(outdir, "input_check.txt")

cat(
  "Input check\n",
  "===========\n",
  file = check_file
)

cat(
  paste0("Ribo-seq count file: ", ribo_file, "\n"),
  paste0("RNA-seq count file: ", rna_file, "\n"),
  paste0("Sample information file: ", data_file, "\n"),
  paste0("Output directory: ", outdir, "\n"),
  paste0("Batch included: ", batch, "\n"),
  paste0("Control condition: ", control, "\n"),
  paste0("Treatment condition: ", treat, "\n"),
  paste0("padj cutoff: ", padj_cutoff, "\n"),
  paste0("log2FC cutoff: ", log2fc_cutoff, "\n"),
  paste0("Ribo-seq sample number: ", ncol(ribo), "\n"),
  paste0("RNA-seq sample number: ", ncol(rna), "\n"),
  paste0("Common gene number: ", length(common_genes), "\n"),
  paste0("Ribo-seq only gene number: ", length(ribo_only_genes), "\n"),
  paste0("RNA-seq only gene number: ", length(rna_only_genes), "\n"),
  paste0("Selected interaction term: ", interaction_term, "\n\n"),
  file = check_file,
  append = TRUE
)

cat(
  "DESeq2 resultsNames:\n",
  paste(res_names, collapse = "\n"),
  "\n\nFinal sample order:\n",
  paste(colnames(merge), collapse = "\n"),
  "\n",
  file = check_file,
  append = TRUE
)

# Choose the interaction term.
# This means:
#   (RIBO treat vs control) - (RNA treat vs control)

res = results(ddsMat, name = interaction_term)

# summary(res)

te_sig = sig_filter(res)
length(which(te_sig))

write.table(
  rownames(res)[which(te_sig)],
  file.path(outdir, "gene_lists", "DTEGs.txt"),
  quote = FALSE,
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE
)

write.table(
  res,
  file.path(outdir, "fold_changes", "deltaTE.txt"),
  quote = FALSE,
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE
)

write_deseq_result(
  res_obj = res,
  outfile = file.path(outdir, "fold_changes", "deltaTE_addsample.txt"),
  control = control,
  treat = treat
)

pdf(file.path(outdir, "Result_figures.pdf"), useDingbats = FALSE)

## Visualisation and interpretation

### DESeq2 object with batch for Ribo-seq

ind = which(coldata$SeqType == "RIBO")
coldata_ribo = coldata[ind, ]
ribo = ribo[, rownames(coldata_ribo), drop = FALSE]

# PCA
if (batch == 1) {
  ddsMat_ribo = DESeqDataSetFromMatrix(
    countData = ribo,
    colData = coldata_ribo,
    design = ~ Condition + Batch
  )
  
  vsd = vst(ddsMat_ribo)
  pcaData = plotPCA(vsd, intgroup = c("Condition_original", "Batch"), returnData = TRUE)
  percentVar = round(100 * attr(pcaData, "percentVar"))
  
  p = ggplot(pcaData, aes(PC1, PC2, color = Condition_original, shape = Batch)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("Ribo-seq PCA") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    theme(aspect.ratio = 1)
  
  print(p)
  
} else if (batch == 0) {
  ddsMat_ribo = DESeqDataSetFromMatrix(
    countData = ribo,
    colData = coldata_ribo,
    design = ~ Condition
  )
  
  vsd = vst(ddsMat_ribo)
  pcaData = plotPCA(vsd, intgroup = "Condition_original", returnData = TRUE)
  percentVar = round(100 * attr(pcaData, "percentVar"))
  
  p = ggplot(pcaData, aes(PC1, PC2, color = Condition_original)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("Ribo-seq PCA") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    theme(aspect.ratio = 1)
  
  print(p)
}

ddsMat_ribo = DESeq(ddsMat_ribo)
ribo_res_names = resultsNames(ddsMat_ribo)
print(ribo_res_names)

ribo_coef = grep("^Condition_treat_vs_control$", ribo_res_names, value = TRUE)

if (length(ribo_coef) != 1) {
  ribo_coef = grep("Condition.*treat.*control", ribo_res_names, value = TRUE)
}

if (length(ribo_coef) != 1) {
  stop(
    paste0(
      "Cannot identify Ribo-seq Condition coefficient.\n",
      "Available resultsNames are:\n",
      paste(ribo_res_names, collapse = "\n")
    ),
    call. = FALSE
  )
}

res_ribo = results(ddsMat_ribo, name = ribo_coef)

#if (requireNamespace("apeglm", quietly = TRUE)) {
#  res_ribo = lfcShrink(ddsMat_ribo, coef = ribo_coef, res = res_ribo, type = "apeglm")
#} else {
#  warning("Package apeglm is not installed. Use unshrunken Ribo-seq log2FoldChange.")
#}

write.table(
  res_ribo,
  file.path(outdir, "fold_changes", "deltaRibo.txt"),
  quote = FALSE,
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE
)

write_deseq_result(
  res_obj = res_ribo,
  outfile = file.path(outdir, "fold_changes", "deltaRibo_addsample.txt"),
  control = control,
  treat = treat
)


### DESeq2 object with batch for RNA-seq

ind = which(coldata$SeqType == "RNA")
coldata_rna = coldata[ind, ]
rna = rna[, rownames(coldata_rna), drop = FALSE]

# PCA
if (batch == 1) {
  ddsMat_rna = DESeqDataSetFromMatrix(
    countData = rna,
    colData = coldata_rna,
    design = ~ Condition + Batch
  )
  
  vsd = vst(ddsMat_rna)
  pcaData = plotPCA(vsd, intgroup = c("Condition_original", "Batch"), returnData = TRUE)
  percentVar = round(100 * attr(pcaData, "percentVar"))
  
  p = ggplot(pcaData, aes(PC1, PC2, color = Condition_original, shape = Batch)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("RNA-seq PCA") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    theme(aspect.ratio = 1)
  
  print(p)
  
} else if (batch == 0) {
  ddsMat_rna = DESeqDataSetFromMatrix(
    countData = rna,
    colData = coldata_rna,
    design = ~ Condition
  )
  
  vsd = vst(ddsMat_rna)
  pcaData = plotPCA(vsd, intgroup = "Condition_original", returnData = TRUE)
  percentVar = round(100 * attr(pcaData, "percentVar"))
  
  p = ggplot(pcaData, aes(PC1, PC2, color = Condition_original)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("RNA-seq PCA") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    theme(aspect.ratio = 1)
  
  print(p)
}

ddsMat_rna = DESeq(ddsMat_rna)
rna_res_names = resultsNames(ddsMat_rna)
print(rna_res_names)

rna_coef = grep("^Condition_treat_vs_control$", rna_res_names, value = TRUE)

if (length(rna_coef) != 1) {
  rna_coef = grep("Condition.*treat.*control", rna_res_names, value = TRUE)
}

if (length(rna_coef) != 1) {
  stop(
    paste0(
      "Cannot identify RNA-seq Condition coefficient.\n",
      "Available resultsNames are:\n",
      paste(rna_res_names, collapse = "\n")
    ),
    call. = FALSE
  )
}

res_rna = results(ddsMat_rna, name = rna_coef)

#if (requireNamespace("apeglm", quietly = TRUE)) {
#  res_rna = lfcShrink(ddsMat_rna, coef = rna_coef, res = res_rna, type = "apeglm")
#} else {
#  warning("Package apeglm is not installed. Use unshrunken RNA-seq log2FoldChange.")
#}

write.table(
  res_rna,
  file.path(outdir, "fold_changes", "deltaRNA.txt"),
  quote = FALSE,
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE
)

write_deseq_result(
  res_obj = res_rna,
  outfile = file.path(outdir, "fold_changes", "deltaRNA_addsample.txt"),
  control = control,
  treat = treat
)



rna_sig = sig_filter(res_rna)

write.table(
  rownames(res_rna)[which(rna_sig)],
  file.path(outdir, "gene_lists", "DEG.txt"),
  quote = FALSE,
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE
)


## Classes of genes

te_sig = sig_filter(res)
ribo_sig = sig_filter(res_ribo)
rna_sig = sig_filter(res_rna)

forwarded = rownames(res)[which(!te_sig & ribo_sig & rna_sig)]
write.table(
  forwarded,
  file.path(outdir, "gene_lists", "forwarded.txt"),
  quote = FALSE,
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE
)

exclusive = rownames(res)[which(te_sig & ribo_sig & !rna_sig)]
write.table(
  exclusive,
  file.path(outdir, "gene_lists", "exclusive.txt"),
  quote = FALSE,
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE
)

both = which(te_sig & ribo_sig & rna_sig)

intensified = rownames(res)[both[which(res[both, "log2FoldChange"] * res_rna[both, "log2FoldChange"] > 0)]]
write.table(
  intensified,
  file.path(outdir, "gene_lists", "intensified.txt"),
  quote = FALSE,
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE
)

buffered = rownames(res)[both[which(res[both, "log2FoldChange"] * res_rna[both, "log2FoldChange"] < 0)]]
buffered = c(
  rownames(res)[which(te_sig & !ribo_sig & rna_sig)],
  buffered
)
buffered = unique(buffered)

write.table(
  buffered,
  file.path(outdir, "gene_lists", "buffered.txt"),
  quote = FALSE,
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE
)

cat(
  "\nGene class summary\n",
  "==================\n",
  paste0("DTEGs: ", length(rownames(res)[which(te_sig)]), "\n"),
  paste0("DEGs: ", length(rownames(res_rna)[which(rna_sig)]), "\n"),
  paste0("Forwarded(!te_sig & ribo_sig & rna_sig): ", length(forwarded), "\n"),
  paste0("Exclusive(te_sig & ribo_sig & !rna_sig): ", length(exclusive), "\n"),
  paste0("Intensified(all sig & te/rna same direct): ", length(intensified), "\n"),
  paste0("Buffered(all sig & te/rna diff direct || te_sig & !ribo_sig & rna_sig ): ", length(buffered), "\n"),
  file = check_file,
  append = TRUE
)


## Scatter plot

max_val = max(abs(c(res_ribo[, "log2FoldChange"], res_rna[, "log2FoldChange"])), na.rm = TRUE)

if (!is.finite(max_val) || max_val == 0) {
  max_val = 1
}

plot(
  y = res_ribo[, "log2FoldChange"],
  x = res_rna[, "log2FoldChange"],
  xlab = "RNA-seq log2 fold change",
  ylab = "Ribo-seq log2 fold change",
  asp = 1,
  pch = 16,
  col = rgb(128 / 255, 128 / 255, 128 / 255, 0.1),
  ylim = c(-max_val, max_val),
  xlim = c(-max_val, max_val),
  cex = 0.4
)

abline(a = 0, b = 1, col = "gray")
abline(h = 0, v = 0, col = "gray")

if (length(forwarded) > 0) {
  points(
    y = res_ribo[forwarded, "log2FoldChange"],
    x = res_rna[forwarded, "log2FoldChange"],
    pch = 16,
    col = rgb(0, 0, 1, 1)
  )
}

if (length(exclusive) > 0) {
  points(
    y = res_ribo[exclusive, "log2FoldChange"],
    x = res_rna[exclusive, "log2FoldChange"],
    pch = 16,
    col = rgb(1, 0, 0, 1)
  )
}

if (length(intensified) > 0) {
  points(
    y = res_ribo[intensified, "log2FoldChange"],
    x = res_rna[intensified, "log2FoldChange"],
    pch = 16,
    col = rgb(1, 0, 1, 1)
  )
}

if (length(buffered) > 0) {
  points(
    y = res_ribo[buffered, "log2FoldChange"],
    x = res_rna[buffered, "log2FoldChange"],
    pch = 16,
    col = rgb(1, 0.5, 0, 1)
  )
}

legend(
  "topleft",
  legend = c("Forwarded", "Exclusive", "Intensified", "Buffered"),
  col = c(
    rgb(0, 0, 1, 1),
    rgb(1, 0, 0, 1),
    rgb(1, 0, 1, 1),
    rgb(1, 0.5, 0, 1)
  ),
  pch = 16,
  bty = "n"
)


### Examples for each class of genes

plot_one_class = function(gene_list, main_title) {
  if (length(gene_list) == 0 || is.na(gene_list[1])) {
    plot.new()
    title(main = main_title)
    text(0.5, 0.5, "No genes in this class")
  } else {
    goi = gene_list[1]
    
    y_u = max(
      res[goi, "log2FoldChange"],
      res_ribo[goi, "log2FoldChange"],
      res_rna[goi, "log2FoldChange"],
      0,
      na.rm = TRUE
    )
    
    y_l = min(
      res[goi, "log2FoldChange"],
      res_ribo[goi, "log2FoldChange"],
      res_rna[goi, "log2FoldChange"],
      0,
      na.rm = TRUE
    )
    
    if (!is.finite(y_u) || !is.finite(y_l) || y_u == y_l) {
      y_u = 1
      y_l = -1
    }
    
    plot(
      c(1, 2),
      c(0, res[goi, "log2FoldChange"]),
      type = "l",
      col = "red",
      xaxt = "n",
      xlab = "Conditions",
      ylim = c(y_l, y_u),
      ylab = "Log2 Fold Change",
      main = paste(main_title, goi, sep = "\n")
    )
    
    lines(c(1, 2), c(0, res_ribo[goi, "log2FoldChange"]), col = "gray")
    lines(c(1, 2), c(0, res_rna[goi, "log2FoldChange"]), col = "blue")
    axis(1, at = c(1, 2), labels = c(control, treat), las = 1)
    
    legend(
      "bottomleft",
      c("RNA", "Ribo", "deltaTE"),
      fill = c("blue", "gray", "red"),
      cex = 1,
      border = NA,
      bty = "n"
    )
  }
}

par(mfrow = c(2, 2))

plot_one_class(forwarded, "Forwarded gene")
plot_one_class(exclusive, "Exclusive gene")
plot_one_class(buffered, "Buffered gene")
plot_one_class(intensified, "Intensified gene")

dev.off()

message("Analysis finished.")
message("Output directory: ", outdir)
message("Main result: ", file.path(outdir, "fold_changes", "deltaTE.txt"))
message("DTEG list: ", file.path(outdir, "gene_lists", "DTEGs.txt"))
message("DEG list: ", file.path(outdir, "gene_lists", "DEG.txt"))