#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(RColorBrewer)
})

# =========================
# 参数设置
# =========================
parser <- ArgumentParser(description = "Expression matrix QC: density plot, sample correlation and PCA")

parser$add_argument("-e", "--expr", required = TRUE,
                    help = "Expression matrix, first column is gene ID, other columns are samples")

parser$add_argument("-g", "--group", required = TRUE,
                    help = "Sample group file, must contain columns: Sample and Group")

parser$add_argument("-o", "--outdir", default = "Expression_QC",
                    help = "Output directory [default: Expression_QC]")

parser$add_argument("--prefix", default = "expr",
                    help = "Output file prefix [default: expr]")

args <- parser$parse_args()

expr_file <- args$expr
group_file <- args$group
outdir <- args$outdir
prefix <- args$prefix

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 读取数据
# =========================
expr <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = "",
  quote = ""
)

group <- read.table(
  group_file,
  header = FALSE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = "",
  quote = ""
)

group <- group[, 1:2, drop = FALSE]
colnames(group) <- c("Group","Sample" )

# 表达矩阵
gene_id <- rownames(expr)
expr_mat <- expr

# 保留原始样品名，防止 Col-3D_1 被自动改成 Col.3D_1
sample_names <- colnames(expr_mat)

# 转为数值矩阵
expr_mat <- as.data.frame(
  lapply(expr_mat, as.numeric),
  check.names = FALSE
)

colnames(expr_mat) <- sample_names
rownames(expr_mat) <- gene_id

# 去除全部为 NA 的基因
expr_mat <- expr_mat[rowSums(is.na(expr_mat)) < ncol(expr_mat), , drop = FALSE]

# NA 替换为 0
expr_mat[is.na(expr_mat)] <- 0

# 检查样品名是否匹配
expr_samples <- colnames(expr_mat)
group_samples <- group$Sample

missing_in_group <- setdiff(expr_samples, group_samples)
missing_in_expr <- setdiff(group_samples, expr_samples)

if (length(missing_in_group) > 0) {
  stop(paste("These samples are in expression matrix but not in group file:",
             paste(missing_in_group, collapse = ", ")))
}

if (length(missing_in_expr) > 0) {
  warning(paste("These samples are in group file but not in expression matrix:",
                paste(missing_in_expr, collapse = ", ")))
}

# 按分组表顺序整理表达矩阵
expr_mat <- expr_mat[, group$Sample, drop = FALSE]

# 过滤低表达基因
expr_mat <- expr_mat[rowSums(expr_mat > 1) > 2, ]

# 更新样品顺序
expr_samples <- colnames(expr_mat)

# log2 转换
log_expr <- log2(expr_mat + 1)

# 输出 log2 表达矩阵
log_expr_out <- cbind(GeneID = rownames(log_expr), log_expr)
write.table(
  log_expr_out,
  file = file.path(outdir, paste0(prefix, ".log2_expression.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================
# 1. 表达量密度分布图
# =========================
density_data <- melt(
  as.matrix(log_expr),
  varnames = c("GeneID", "Sample"),
  value.name = "Expression"
)

density_data <- merge(density_data, group, by = "Sample", all.x = TRUE)

p_density <- ggplot(density_data, aes(x = Expression, color = Sample)) +
  geom_density(linewidth = 0.5) +
  theme_bw(base_size = 13) +
  labs(
    title = "Expression density distribution",
    x = "log2(Expression + 1)",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

ggsave(
  filename = file.path(outdir, paste0(prefix, ".density.pdf")),
  plot = p_density,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(outdir, paste0(prefix, ".density.png")),
  plot = p_density,
  width = 8,
  height = 6,
  dpi = 300
)

# =========================
# 1.1 表达量箱线图
# =========================
p_boxplot <- ggplot(density_data, aes(x = Sample, y = Expression, fill = Group)) +
  geom_boxplot(
    outlier.size = 0.3,
    linewidth = 0.4,
    width = 0.5
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Expression boxplot",
    x = "Sample",
    y = "log2(Expression + 1)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(outdir, paste0(prefix, ".boxplot.pdf")),
  plot = p_boxplot,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(outdir, paste0(prefix, ".boxplot.png")),
  plot = p_boxplot,
  width = 8,
  height = 6,
  dpi = 300
)

# =========================
# 2. 样品间相关性分析
# =========================
cor_mat <- cor(log_expr, method = "pearson", use = "pairwise.complete.obs")

write.table(
  round(cor_mat, 4),
  file = file.path(outdir, paste0(prefix, ".sample_correlation.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

annotation_col <- data.frame(Group = group$Group)
rownames(annotation_col) <- group$Sample

group_levels <- unique(group$Group)
ann_colors <- NULL

if (length(group_levels) <= 12) {
  cols <- brewer.pal(max(3, min(12, length(group_levels))), "Set3")[seq_along(group_levels)]
  names(cols) <- group_levels
  ann_colors <- list(Group = cols)
}

pdf(file.path(outdir, paste0(prefix, ".sample_correlation_heatmap.pdf")),
    width = 7,
    height = 6)

pheatmap(
  cor_mat,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = ann_colors,
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows = F,
  cluster_cols = F,
  #clustering_distance_rows = "euclidean",
  #clustering_distance_cols = "euclidean",
  #clustering_method = "complete",
  main = "Sample correlation",
  fontsize = 10
)

dev.off()

png(file.path(outdir, paste0(prefix, ".sample_correlation_heatmap.png")),
    width = 2100,
    height = 1800,
    res = 300)

pheatmap(
  cor_mat,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = ann_colors,
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows = F,
  cluster_cols = F,
  #clustering_distance_rows = "euclidean",
  #clustering_distance_cols = "euclidean",
  #clustering_method = "complete",
  main = "Sample correlation",
  fontsize = 10
)

dev.off()

# =========================
# 3. PCA 分析
# =========================

# 去除方差为 0 的基因
gene_var <- apply(log_expr, 1, var)
pca_mat <- log_expr[gene_var > 0, , drop = FALSE]

pca_res <- prcomp(t(pca_mat), center = TRUE, scale. = FALSE)

pca_data <- as.data.frame(pca_res$x)
pca_data$Sample <- rownames(pca_data)
pca_data <- merge(pca_data, group, by = "Sample", all.x = TRUE)

var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

write.table(
  pca_data,
  file = file.path(outdir, paste0(prefix, ".PCA_coordinates.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
#  geom_text(vjust = -0.8, size = 3.5, show.legend = FALSE) +
  theme_bw(base_size = 13) +
  labs(
    title = "PCA plot",
    x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

ggsave(
  filename = file.path(outdir, paste0(prefix, ".PCA.pdf")),
  plot = p_pca,
  width = 7,
  height = 6
)

ggsave(
  filename = file.path(outdir, paste0(prefix, ".PCA.png")),
  plot = p_pca,
  width = 7,
  height = 6,
  dpi = 300
)

# =========================
# 输出完成信息
# =========================
cat("Expression QC finished.\n")
cat("Output directory:", outdir, "\n")