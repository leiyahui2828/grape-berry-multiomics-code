#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# =========================
# 1. 参数设置
# =========================

option_list <- list(
  make_option(c("--ribo"), type = "character",
              help = "Ribo-seq expression matrix, first column is gene ID"),
  make_option(c("--rna"), type = "character",
              help = "RNA-seq expression matrix, first column is gene ID"),
  make_option(c("--sample_info"), type = "character",
              help = "Sample information file with header: SampleID Condition SeqType ..."),
  make_option(c("--outdir"), type = "character",
              default = "ribo_rna_correlation_result",
              help = "Output directory [default: %default]"),
  make_option(c("--method"), type = "character",
              default = "pearson",
              help = "Correlation method for heatmap: pearson, spearman or kendall [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$ribo) || is.null(opt$rna) || is.null(opt$sample_info)) {
  stop("Please provide --ribo, --rna and --sample_info")
}

if (!opt$method %in% c("pearson", "spearman", "kendall")) {
  stop("--method must be one of: pearson, spearman, kendall")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

mean_scatter_dir <- file.path(opt$outdir, "condition_mean_scatter_plots")
dir.create(mean_scatter_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 2. 读取表达矩阵
# =========================

read_expr_matrix <- function(file) {
  dat <- read.table(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  
  if (ncol(dat) < 2) {
    stop(paste("Expression matrix should contain at least two columns:", file))
  }
  
  gene_id <- dat[[1]]
  expr <- dat[, -1, drop = FALSE]
  
  sample_names <- colnames(expr)
  
  expr <- as.data.frame(
    lapply(expr, function(x) as.numeric(as.character(x))),
    check.names = FALSE
  )
  
  colnames(expr) <- sample_names
  rownames(expr) <- gene_id
  
  keep <- !is.na(rownames(expr)) & rownames(expr) != ""
  expr <- expr[keep, , drop = FALSE]
  
  if (any(duplicated(rownames(expr)))) {
    message("Detected duplicated gene IDs in ", file, ", keeping the first occurrence.")
    expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  }
  
  return(expr)
}

ribo_raw <- read_expr_matrix(opt$ribo)
rna_raw  <- read_expr_matrix(opt$rna)

# =========================
# 3. 读取 sample_info.txt
# =========================

sample_info <- read.table(
  opt$sample_info,
  header = TRUE,
  sep = "",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  blank.lines.skip = TRUE
)

required_cols <- c("SampleID", "Condition", "SeqType")
missing_cols <- setdiff(required_cols, colnames(sample_info))

if (length(missing_cols) > 0) {
  stop(
    "sample_info file must contain columns: SampleID, Condition, SeqType. Missing: ",
    paste(missing_cols, collapse = ", ")
  )
}

sample_info <- sample_info[, required_cols]
sample_info$SeqType <- toupper(sample_info$SeqType)

sample_info <- sample_info[sample_info$SeqType %in% c("RIBO", "RNA"), , drop = FALSE]

if (nrow(sample_info) == 0) {
  stop("No valid samples found. SeqType should be RIBO or RNA.")
}

ribo_info <- sample_info[sample_info$SeqType == "RIBO", , drop = FALSE]
rna_info  <- sample_info[sample_info$SeqType == "RNA",  , drop = FALSE]

if (nrow(ribo_info) == 0) {
  stop("No RIBO samples found in sample_info.")
}

if (nrow(rna_info) == 0) {
  stop("No RNA samples found in sample_info.")
}

missing_ribo <- setdiff(ribo_info$SampleID, colnames(ribo_raw))
missing_rna  <- setdiff(rna_info$SampleID, colnames(rna_raw))

if (length(missing_ribo) > 0) {
  stop(
    "The following RIBO samples are missing in Ribo-seq matrix: ",
    paste(missing_ribo, collapse = ", ")
  )
}

if (length(missing_rna) > 0) {
  stop(
    "The following RNA samples are missing in RNA-seq matrix: ",
    paste(missing_rna, collapse = ", ")
  )
}

# 按 sample_info 中的顺序排序
ribo_raw <- ribo_raw[, ribo_info$SampleID, drop = FALSE]
rna_raw  <- rna_raw[, rna_info$SampleID, drop = FALSE]

# 不添加 RIBO/RNA 前缀时，要求所有样品名唯一
combined_sample_names <- c(colnames(ribo_raw), colnames(rna_raw))

if (any(duplicated(combined_sample_names))) {
  dup_samples <- unique(combined_sample_names[duplicated(combined_sample_names)])
  stop(
    "Duplicated sample names found between RIBO and RNA matrices: ",
    paste(dup_samples, collapse = ", "),
    "\nBecause sample names are displayed without RIBO/RNA prefix, ",
    "RIBO and RNA sample names must be unique."
  )
}

# =========================
# 4. 取共有基因
# =========================

common_genes <- intersect(rownames(ribo_raw), rownames(rna_raw))

if (length(common_genes) == 0) {
  stop("No common genes found between Ribo-seq and RNA-seq matrices.")
}

message("Number of common genes: ", length(common_genes))

ribo <- ribo_raw[common_genes, , drop = FALSE]
rna  <- rna_raw[common_genes, , drop = FALSE]

if (any(ribo < 0, na.rm = TRUE) || any(rna < 0, na.rm = TRUE)) {
  stop("Expression matrix contains negative values. log2(expression + 1) cannot be applied.")
}

# =========================
# 5. 所有样品相关性热图和矩阵
#    单样品层面：先 log2(expression + 1)，再计算相关性
#    样品名不额外添加 RIBO/RNA 前缀
# =========================

ribo_log <- log2(ribo + 1)
rna_log  <- log2(rna + 1)

combined_expr <- cbind(ribo_log, rna_log)

cor_mat <- cor(
  combined_expr,
  method = opt$method,
  use = "pairwise.complete.obs"
)

write.table(
  round(cor_mat, 4),
  file = file.path(opt$outdir, "01.all_sample_correlation_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

annotation_col <- data.frame(
  SeqType = c(ribo_info$SeqType, rna_info$SeqType),
  Condition = c(ribo_info$Condition, rna_info$Condition),
  stringsAsFactors = FALSE
)

rownames(annotation_col) <- colnames(combined_expr)

condition_levels <- unique(sample_info$Condition)

if (length(condition_levels) <= 8) {
  condition_cols <- brewer.pal(max(3, length(condition_levels)), "Set2")[seq_along(condition_levels)]
} else {
  condition_cols <- colorRampPalette(brewer.pal(8, "Set2"))(length(condition_levels))
}
names(condition_cols) <- condition_levels

annotation_colors <- list(
  SeqType = c("RIBO" = "#E64B35", "RNA" = "#4DBBD5"),
  Condition = condition_cols
)

pdf(
  file.path(opt$outdir, "02.all_sample_correlation_heatmap.pdf"),
  width = max(10, ncol(combined_expr) * 0.5+1),
  height = max(8, ncol(combined_expr) * 0.5)
)

pheatmap(
  cor_mat,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  #breaks = seq(-1, 1, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = annotation_colors,
  #fontsize_number = 10, 
  fontsize = 12,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = paste0("All sample correlation heatmap (", opt$method, ")")
)

dev.off()

# =========================
# 6. 计算 Condition 组内均值表达量
#    组均值层面：先在原始表达量尺度取均值，再 log2(mean + 1)
# =========================

calc_condition_mean <- function(expr_raw, info) {
  conditions <- unique(info$Condition)
  
  mean_mat <- sapply(conditions, function(condition_name) {
    samples <- info$SampleID[info$Condition == condition_name]
    samples <- intersect(samples, colnames(expr_raw))
    
    if (length(samples) == 0) {
      return(rep(NA_real_, nrow(expr_raw)))
    } else if (length(samples) == 1) {
      return(expr_raw[, samples])
    } else {
      return(rowMeans(expr_raw[, samples, drop = FALSE], na.rm = TRUE))
    }
  })
  
  mean_mat <- as.data.frame(mean_mat, check.names = FALSE)
  rownames(mean_mat) <- rownames(expr_raw)
  
  return(mean_mat)
}

ribo_mean_raw <- calc_condition_mean(ribo, ribo_info)
rna_mean_raw  <- calc_condition_mean(rna, rna_info)

ribo_mean <- log2(ribo_mean_raw + 1)
rna_mean  <- log2(rna_mean_raw + 1)

shared_conditions <- intersect(colnames(ribo_mean), colnames(rna_mean))

if (length(shared_conditions) == 0) {
  warning("No shared Condition found between RIBO and RNA samples. Skip condition-mean analysis.")
} else {
  
  ribo_mean <- ribo_mean[, shared_conditions, drop = FALSE]
  rna_mean  <- rna_mean[, shared_conditions, drop = FALSE]
  
  # =========================
  # 7. Condition 均值相关性矩阵和热图
  # =========================
  
  ribo_mean_for_cor <- ribo_mean
  rna_mean_for_cor  <- rna_mean
  
  # 均值矩阵中同一个 Condition 同时有 RIBO 和 RNA，必须区分列名
  colnames(ribo_mean_for_cor) <- paste0("RIBO_", colnames(ribo_mean_for_cor))
  colnames(rna_mean_for_cor)  <- paste0("RNA_",  colnames(rna_mean_for_cor))
  
  condition_mean_expr <- cbind(ribo_mean_for_cor, rna_mean_for_cor)
  
  condition_mean_cor_mat <- cor(
    condition_mean_expr,
    method = opt$method,
    use = "pairwise.complete.obs"
  )
  
  write.table(
    round(condition_mean_cor_mat, 4),
    file = file.path(opt$outdir, "03.condition_mean_correlation_matrix.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  condition_mean_anno <- data.frame(
    SeqType = c(
      rep("RIBO", ncol(ribo_mean_for_cor)),
      rep("RNA", ncol(rna_mean_for_cor))
    ),
    Condition = c(shared_conditions, shared_conditions),
    stringsAsFactors = FALSE
  )
  
  rownames(condition_mean_anno) <- colnames(condition_mean_expr)
  
  mean_annotation_colors <- list(
    SeqType = c("RIBO" = "#E64B35", "RNA" = "#4DBBD5"),
    Condition = condition_cols[names(condition_cols) %in% shared_conditions]
  )
  
  pdf(
    file.path(opt$outdir, "04.condition_mean_correlation_heatmap.pdf"),
    width = max(10, ncol(condition_mean_expr) * 0.6 + 1),
    height = max(8, ncol(condition_mean_expr) * 0.6)
  )
  
  pheatmap(
    condition_mean_cor_mat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    #breaks = seq(-1, 1, length.out = 101),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = condition_mean_anno,
    annotation_row = condition_mean_anno,
    annotation_colors = mean_annotation_colors,
    fontsize = 15,
    #fontsize_number = 15, 
    display_numbers = TRUE,
    number_format = "%.2f",
    main = paste0("Condition mean correlation heatmap (", opt$method, ")")
  )
  
  dev.off()
  
  # =========================
  # 8. 每个 Condition 的 RIBO_mean vs RNA_mean 相关性统计和散点图
  # =========================
  

  mean_cor_summary <- data.frame()
  
  for (condition_name in shared_conditions) {
    
    x <- rna_mean[, condition_name]
    y <- ribo_mean[, condition_name]
    
    plot_dat <- data.frame(
      gene_id = common_genes,
      Condition = condition_name,
      RNA = x,
      RIBO = y,
      stringsAsFactors = FALSE
    )
    
    plot_dat <- plot_dat[complete.cases(plot_dat[, c("RNA", "RIBO")]), ]
    
    if (nrow(plot_dat) < 3) {
      warning(
        "Less than 3 valid genes for condition mean: ",
        condition_name,
        ". Skip this condition."
      )
      next
    }
    
    pearson_cor <- cor(
      plot_dat$RNA,
      plot_dat$RIBO,
      method = "pearson",
      use = "complete.obs"
    )
    
    spearman_cor <- cor(
      plot_dat$RNA,
      plot_dat$RIBO,
      method = "spearman",
      use = "complete.obs"
    )
    
    kendall_cor <- cor(
      plot_dat$RNA,
      plot_dat$RIBO,
      method = "kendall",
      use = "complete.obs"
    )
    
    mean_cor_summary <- rbind(
      mean_cor_summary,
      data.frame(
        Condition = condition_name,
        common_gene_number = nrow(plot_dat),
        RIBO_sample_number = length(ribo_info$SampleID[ribo_info$Condition == condition_name]),
        RNA_sample_number = length(rna_info$SampleID[rna_info$Condition == condition_name]),
        pearson_cor = pearson_cor,
        spearman_cor = spearman_cor,
        kendall_cor = kendall_cor,
        stringsAsFactors = FALSE
      )
    )
    
    p <- ggplot(plot_dat, aes(x = RNA, y = RIBO)) +
      geom_point(size = 0.8, alpha = 0.5, color = "grey30") +
      geom_smooth(method = "lm", se = FALSE, color = "#E64B35", linewidth = 0.7) +
      theme_bw(base_size = 12) +
      labs(
        #title = paste0("Condition mean: ", condition_name),
        title =  condition_name,
        subtitle = paste0(
          "Pearson r = ", sprintf("%.3f", pearson_cor),
          " | Spearman rho = ", sprintf("%.3f", spearman_cor)
        ),
        x = "RNA-seq log2(expression + 1)",
        y = "Ribo-seq log2(expression + 1)"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid = element_blank()
      )
    

    ggsave(
      filename = file.path(
        mean_scatter_dir,
        paste0(condition_name, ".condition_mean.ribo_vs_rna.scatter.pdf")
      ),
      plot = p,
      width = 6.5,
      height = 5.5
    )
  }
  
  if (nrow(mean_cor_summary) > 0) {
    mean_cor_summary$pearson_cor  <- round(mean_cor_summary$pearson_cor, 4)
    mean_cor_summary$spearman_cor <- round(mean_cor_summary$spearman_cor, 4)
    mean_cor_summary$kendall_cor  <- round(mean_cor_summary$kendall_cor, 4)
  }
  
  write.table(
    mean_cor_summary,
    file = file.path(opt$outdir, "05.condition_mean_ribo_vs_rna_correlation_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
 }

# =========================
# 9. 输出运行信息
# =========================

message("Analysis finished.")
message("Output directory: ", opt$outdir)
message("Main output files:")
message("  01.all_sample_correlation_matrix.tsv")
message("  02.all_sample_correlation_heatmap.pdf")
message("  03.condition_mean_correlation_matrix.tsv")
message("  04.condition_mean_correlation_heatmap.pdf")
message("  05.condition_mean_ribo_vs_rna_correlation_summary.tsv")
message("  condition_mean_scatter_plots/")