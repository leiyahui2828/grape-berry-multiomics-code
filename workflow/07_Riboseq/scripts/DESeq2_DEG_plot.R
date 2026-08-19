#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(EnhancedVolcano)
  library(ggplot2)
  library(pheatmap)
})

# =========================
# 参数设置
# =========================
parser <- ArgumentParser(
  description = "Extract DEGs and plot volcano/expression heatmap from DESeq2 results"
)

parser$add_argument("-i", "--input", required = TRUE,
                    help = "DESeq2 result file")

parser$add_argument("-e", "--expr", required = TRUE,
                    help = "Expression matrix, first column is GeneID, other columns are samples")

parser$add_argument("-o", "--outdir", default = "DESeq2_DEG_plot",
                    help = "Output directory [default: DESeq2_DEG_plot]")

parser$add_argument("--padj", type = "double", default = 0.05,
                    help = "Adjusted p-value cutoff [default: 0.05]")

parser$add_argument("--log2fc", type = "double", default = 1,
                    help = "Absolute log2 fold change cutoff [default: 1]")

parser$add_argument("--top", type = "integer", default = 50,
                    help = "Top significant genes per comparison for expression heatmap [default: 50]")

parser$add_argument("--label_top", type = "integer", default = 20,
                    help = "Top significant genes to label in volcano plot [default: 20]")

args <- parser$parse_args()

input_file <- args$input
expr_file <- args$expr
outdir <- args$outdir
padj_cutoff <- args$padj
log2fc_cutoff <- args$log2fc
top_n <- args$top
label_top <- args$label_top

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 读取 DESeq2 结果
# =========================
deg <- read.table(
  input_file,
  header = TRUE,
  sep = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = "",
  quote = ""
)

required_cols <- c(
  "GeneID", "sampleA", "sampleB", "baseMeanA", "baseMeanB",
  "log2FoldChange", "pvalue", "padj"
)

missing_cols <- setdiff(required_cols, colnames(deg))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

deg$GeneID <- as.character(deg$GeneID)
deg$sampleA <- as.character(deg$sampleA)
deg$sampleB <- as.character(deg$sampleB)

deg$baseMeanA <- as.numeric(deg$baseMeanA)
deg$baseMeanB <- as.numeric(deg$baseMeanB)
deg$log2FoldChange <- as.numeric(deg$log2FoldChange)
deg$pvalue <- as.numeric(deg$pvalue)
deg$padj <- as.numeric(deg$padj)

# =========================
# 方向检查：仅对显著差异基因
# =========================
check_idx <- !is.na(deg$baseMeanA) &
  !is.na(deg$baseMeanB) &
  !is.na(deg$log2FoldChange) &
  !is.na(deg$padj) &
  deg$padj < padj_cutoff &
  abs(deg$log2FoldChange) >= log2fc_cutoff &
  deg$baseMeanA != deg$baseMeanB &
  deg$log2FoldChange != 0

bad_idx <- check_idx &
  sign(deg$baseMeanA - deg$baseMeanB) != sign(deg$log2FoldChange)

if (any(bad_idx)) {
  bad_file <- file.path(outdir, "log2FoldChange_direction_check.failed.tsv")
  write.table(deg[bad_idx, ], file = bad_file, sep = "\t", quote = FALSE, row.names = FALSE)
  stop(paste0("Direction check failed in significant DEGs.\n",
              "Failed records written to: ", bad_file))
}

# =========================
# 读取表达矩阵
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

gene_id <- rownames(expr)
sample_names <- colnames(expr)

expr_mat <- as.data.frame(lapply(expr, as.numeric), check.names = FALSE)
colnames(expr_mat) <- sample_names
rownames(expr_mat) <- gene_id
expr_mat[is.na(expr_mat)] <- 0
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]
log_expr_mat <- log2(expr_mat + 1)

# =========================
# 生成比较名称
# =========================
deg <- deg[!is.na(deg$sampleA) & !is.na(deg$sampleB) &
             deg$sampleA != "" & deg$sampleB != "", , drop = FALSE]

deg$Comparison <- paste(deg$sampleA, "vs", deg$sampleB, sep = "_")
comparisons <- unique(deg$Comparison)

write.table(deg, file = file.path(outdir, "all.DESeq2.results.with_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_list <- list()
heatmap_gene_list <- c()

for (cmp in comparisons) {
  message("Processing comparison: ", cmp)
  dat <- deg[deg$Comparison == cmp, , drop = FALSE]
  dat <- dat[!is.na(dat$log2FoldChange) & !is.na(dat$padj), , drop = FALSE]
  if (nrow(dat) == 0) next
  
  # 差异状态
  dat$Regulation <- "Not_significant"
  dat$Regulation[dat$padj < padj_cutoff & dat$log2FoldChange >= log2fc_cutoff] <- "Up"
  dat$Regulation[dat$padj < padj_cutoff & dat$log2FoldChange <= -log2fc_cutoff] <- "Down"
  
  sig <- dat[dat$Regulation %in% c("Up","Down"), , drop=FALSE]
  up <- dat[dat$Regulation == "Up", , drop=FALSE]
  down <- dat[dat$Regulation == "Down", , drop=FALSE]
  
  cmp_prefix <- gsub("[^A-Za-z0-9_.-]", "_", cmp)
  
  write.table(dat, file=file.path(outdir, paste0(cmp_prefix,".all.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(sig, file=file.path(outdir, paste0(cmp_prefix,".DEG.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(up, file=file.path(outdir, paste0(cmp_prefix,".DEG.up.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(down, file=file.path(outdir, paste0(cmp_prefix,".DEG.down.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
  
  summary_list[[cmp]] <- data.frame(
    Comparison = cmp,
    Total_genes = nrow(dat),
    DEG_total = nrow(sig),
    Up = nrow(up),
    Down = nrow(down),
    padj_cutoff = padj_cutoff,
    log2fc_cutoff = log2fc_cutoff,
    stringsAsFactors = FALSE
  )
  
  # 收集热图基因
  if (nrow(sig) > 0) {
    sig_ordered <- sig[order(sig$padj), , drop = FALSE]
    heatmap_gene_list <- unique(c(heatmap_gene_list, head(sig_ordered$GeneID, top_n)))
  }
  
  # 火山图
  plot_dat <- dat[!is.na(dat$padj) & !is.na(dat$log2FoldChange), , drop=FALSE]
  plot_dat$padj[plot_dat$padj==0] <- min(plot_dat$padj[plot_dat$padj>0], na.rm=TRUE)
  top_label_genes <- if (nrow(sig)>0 && label_top>0) head(sig[order(sig$padj), "GeneID"], label_top) else character(0)
  volcano_labels <- ifelse(plot_dat$GeneID %in% top_label_genes, plot_dat$GeneID, "")
  
  p_volcano <- EnhancedVolcano(
    plot_dat,
    lab = volcano_labels,
    #lab = plot_dat$GeneID ,
    #selectLab = top_label_genes,
    x = "log2FoldChange",
    y = "padj",
    title = cmp,
    subtitle = paste0("padj < ", padj_cutoff, ", |log2FC| >= ", log2fc_cutoff),
    caption = paste0("Up: ", nrow(up),"  Down: ", nrow(down)),
    pCutoff = padj_cutoff,
    FCcutoff = log2fc_cutoff,
    xlab = bquote(~Log[2]~ "fold change"),
    ylab = bquote(~-Log[10]~ "adjusted P value"),
    pointSize = 1,
    labSize = 3,
    col = c("grey70","grey70","#377EB8","#E41A1C"),
    colAlpha = 0.8,
    legendPosition = "none",
    #legendLabels = c("NS", expression(Log[2]~FC), "padj", expression(padj~and~Log[2]~FC)),
    #legendPosition = "right",
    #legendLabSize = 10,
    #legendIconSize = 3,
    drawConnectors = TRUE,
    widthConnectors = 0.4,
    colConnectors = "grey40",
    max.overlaps = Inf
  )
  
  ggsave(file.path(outdir, paste0(cmp_prefix,".volcano.pdf")), plot=p_volcano, width=8, height=7)
  ggsave(file.path(outdir, paste0(cmp_prefix,".volcano.png")), plot=p_volcano, width=8, height=7, dpi=300)
}

# DEG 统计表
if (length(summary_list) > 0) {
  summary_tab <- do.call(rbind, summary_list)
  write.table(summary_tab, file=file.path(outdir,"DEG.summary.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  
  # =========================
  # 每组比较上调/下调 DEG 数量柱状图
  # =========================
  deg_count_plot_data <- rbind(
    data.frame(
      Comparison = summary_tab$Comparison,
      Regulation = "Up",
      Count = summary_tab$Up,
      stringsAsFactors = FALSE
    ),
    data.frame(
      Comparison = summary_tab$Comparison,
      Regulation = "Down",
      Count = -summary_tab$Down,
      stringsAsFactors = FALSE
    )
  )
  
  deg_count_plot_data$Comparison <- factor(
    deg_count_plot_data$Comparison,
    levels = summary_tab$Comparison
  )
  
  deg_count_plot_data$Regulation <- factor(
    deg_count_plot_data$Regulation,
    levels = c("Up", "Down")
  )
  
  p_deg_count <- ggplot(
    deg_count_plot_data,
    aes(x = Comparison, y = Count, fill = Regulation)
  ) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_y_continuous(labels = abs) +
    scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
    labs(
      title = "DEG number per comparison",
      x = "Comparison",
      y = "Number of DEGs"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    file.path(outdir, "DEG.count.barplot.pdf"),
    plot = p_deg_count,
    width = max(7, nrow(summary_tab) * 0.6),
    height = 5
  )
  
  ggsave(
    file.path(outdir, "DEG.count.barplot.png"),
    plot = p_deg_count,
    width = max(7, nrow(summary_tab) * 0.6),
    height = 5,
    dpi = 300
  )
}

# =========================
# 显著基因表达量热图
# =========================
if (length(heatmap_gene_list) > 1) {
  heatmap_genes <- intersect(unique(heatmap_gene_list), rownames(log_expr_mat))
  if (length(heatmap_genes) > 1) {
    heatmap_mat <- log_expr_mat[heatmap_genes,,drop=FALSE]
    
    write.table(cbind(GeneID=rownames(heatmap_mat), heatmap_mat),
                file=file.path(outdir,"DEG.expression.heatmap_matrix.log2.tsv"),
                sep="\t", quote=FALSE, row.names=FALSE)
    
    heatmap_colors <- colorRampPalette(c("#377EB8","white","#E41A1C"))(100)
    
    pdf(file.path(outdir,"DEG.expression.heatmap.pdf"),
        width=max(7,ncol(heatmap_mat)*0.45),
        height=max(6,nrow(heatmap_mat)*0.15))
    pheatmap(heatmap_mat, scale="row", color=heatmap_colors,
             cluster_rows=TRUE, cluster_cols=TRUE,
             border_color=NA, fontsize_row=ifelse(nrow(heatmap_mat)>80,5,7),
             fontsize_col=10, main="DEG expression heatmap")
    dev.off()
    
    png(file.path(outdir,"DEG.expression.heatmap.png"),
        width=max(2100,ncol(heatmap_mat)*180),
        height=max(1800,nrow(heatmap_mat)*45), res=300)
    pheatmap(heatmap_mat, scale="row", color=heatmap_colors,
             cluster_rows=TRUE, cluster_cols=TRUE,
             border_color=NA, fontsize_row=ifelse(nrow(heatmap_mat)>80,5,7),
             fontsize_col=10, main="DEG top: Z-score of log2(exp + 1)")
    dev.off()
  } else {
    warning("Not enough DEG genes found in expression matrix for heatmap.")
  }
} else {
  warning("No enough significant DEGs for expression heatmap.")
}

message("Finished.")
message("Output directory: ", outdir)
