#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
})

# ============================================================
# 1. Parameters
# ============================================================

option_list <- list(
  make_option(
    c("--ribo_expr"),
    type = "character",
    help = "Ribo-seq expression matrix"
  ),
  make_option(
    c("--rna_expr"),
    type = "character",
    help = "RNA-seq expression matrix"
  ),
  make_option(
    c("--sample_info"),
    type = "character",
    help = "Sample information file with columns: SampleID, Condition, SeqType"
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = "TE_out",
    help = "Output directory [default: %default]"
  ),
  make_option(
    c("--mean_cutoff"),
    type = "double",
    default = 1,
    help = "RNA group mean cutoff [default: %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$ribo_expr) || is.null(opt$rna_expr) || is.null(opt$sample_info)) {
  stop("Please provide --ribo_expr, --rna_expr and --sample_info")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. Read expression matrix
#    First column is treated as GeneID
# ============================================================

read_expr_matrix <- function(file) {
  
  dat <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = ""
  )
  
  if (ncol(dat) < 2) {
    stop("Expression matrix should contain at least two columns: ", file)
  }
  
  colnames(dat)[1] <- "GeneID"
  dat$GeneID <- as.character(dat$GeneID)
  
  if (anyDuplicated(dat$GeneID) > 0) {
    dup_gene <- dat$GeneID[duplicated(dat$GeneID)][1]
    stop("Duplicated GeneID found in ", file, ": ", dup_gene)
  }
  
  expr <- dat %>%
    column_to_rownames("GeneID")
  
  expr <- as.data.frame(expr, check.names = FALSE)
  expr[] <- lapply(expr, as.numeric)
  
  return(expr)
}


# ============================================================
# 3. Read sample_info
#    Required columns: SampleID, Condition, SeqType
# ============================================================

read_sample_info <- function(file) {
  
  info <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = ""
  )
  
  required_cols <- c("SampleID", "Condition", "SeqType")
  missing_cols <- setdiff(required_cols, colnames(info))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in sample_info: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  info <- info %>%
    select(SampleID, Condition, SeqType) %>%
    mutate(
      SampleID = as.character(SampleID),
      Condition = as.character(Condition),
      SeqType = toupper(as.character(SeqType))
    )
  
  invalid_seqtype <- setdiff(unique(info$SeqType), c("RIBO", "RNA"))
  
  if (length(invalid_seqtype) > 0) {
    stop(
      "Invalid SeqType found in sample_info: ",
      paste(invalid_seqtype, collapse = ", "),
      ". SeqType should be RIBO or RNA."
    )
  }
  
  if (anyDuplicated(info$SampleID) > 0) {
    dup_sample <- info$SampleID[duplicated(info$SampleID)][1]
    stop("Duplicated SampleID found in sample_info: ", dup_sample)
  }
  
  return(info)
}


# ============================================================
# 4. Calculate group mean expression by Condition
# ============================================================

calculate_group_mean <- function(expr, sample_info, seqtype, expr_name = "expr") {
  
  sample_sub <- sample_info %>%
    filter(SeqType == seqtype)
  
  if (nrow(sample_sub) == 0) {
    stop("No samples with SeqType = ", seqtype, " in sample_info.")
  }
  
  missing_samples <- setdiff(sample_sub$SampleID, colnames(expr))
  
  if (length(missing_samples) > 0) {
    stop(
      "Samples in sample_info not found in ",
      expr_name,
      " matrix: ",
      paste(missing_samples, collapse = ", ")
    )
  }
  
  conditions <- unique(sample_sub$Condition)
  
  group_mean <- map_dfc(
    conditions,
    function(cond) {
      
      samples <- sample_sub %>%
        filter(Condition == cond) %>%
        pull(SampleID)
      
      if (length(samples) == 1) {
        tibble(!!cond := expr[, samples])
      } else {
        tibble(!!cond := rowMeans(expr[, samples, drop = FALSE], na.rm = TRUE))
      }
    }
  )
  
  group_mean <- as.data.frame(group_mean, check.names = FALSE)
  rownames(group_mean) <- rownames(expr)
  
  return(group_mean)
}


# ============================================================
# 5. Read input data
# ============================================================

ribo_expr <- read_expr_matrix(opt$ribo_expr)
rna_expr <- read_expr_matrix(opt$rna_expr)
sample_info <- read_sample_info(opt$sample_info)

common_genes <- intersect(rownames(ribo_expr), rownames(rna_expr))

if (length(common_genes) == 0) {
  stop("No common genes found between Ribo-seq and RNA-seq expression matrices.")
}

ribo_expr <- ribo_expr[common_genes, , drop = FALSE]
rna_expr <- rna_expr[common_genes, , drop = FALSE]


# ============================================================
# 6. Check common conditions
# ============================================================

ribo_conditions <- sample_info %>%
  filter(SeqType == "RIBO") %>%
  distinct(Condition) %>%
  pull(Condition)

rna_conditions <- sample_info %>%
  filter(SeqType == "RNA") %>%
  distinct(Condition) %>%
  pull(Condition)

common_conditions <- intersect(ribo_conditions, rna_conditions)

if (length(common_conditions) == 0) {
  stop("No common Condition found between RIBO and RNA samples.")
}

condition_status <- full_join(
  tibble(Condition = ribo_conditions, in_ribo = TRUE),
  tibble(Condition = rna_conditions, in_rna = TRUE),
  by = "Condition"
) %>%
  mutate(
    in_ribo = replace_na(in_ribo, FALSE),
    in_rna = replace_na(in_rna, FALSE),
    condition_status = case_when(
      in_ribo & in_rna ~ "Both",
      in_ribo & !in_rna ~ "RIBO_only",
      !in_ribo & in_rna ~ "RNA_only",
      TRUE ~ "Unknown"
    )
  )

write.table(
  condition_status,
  file = file.path(opt$outdir, "condition_status.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sample_info_valid <- sample_info %>%
  filter(Condition %in% common_conditions)


# ============================================================
# 7. Calculate group mean expression
# ============================================================

ribo_group_mean <- calculate_group_mean(
  expr = ribo_expr,
  sample_info = sample_info_valid,
  seqtype = "RIBO",
  expr_name = "Ribo-seq"
)

rna_group_mean <- calculate_group_mean(
  expr = rna_expr,
  sample_info = sample_info_valid,
  seqtype = "RNA",
  expr_name = "RNA-seq"
)

ribo_group_mean <- ribo_group_mean[, common_conditions, drop = FALSE]
rna_group_mean <- rna_group_mean[, common_conditions, drop = FALSE]


# ============================================================
# 8. Calculate TE
#    TE = Ribo group mean / RNA group mean
#    Genes with RNA/Ribo group mean < cutoff are set to NA
# ============================================================

te_matrix <- ribo_group_mean / rna_group_mean

te_matrix[rna_group_mean < opt$mean_cutoff] <- NA
#te_matrix[ribo_group_mean < opt$mean_cutoff] <- NA

te_matrix[is.infinite(as.matrix(te_matrix))] <- NA

log2te_matrix <- log2(te_matrix)
log2te_matrix[is.infinite(as.matrix(log2te_matrix))] <- NA


# ============================================================
# 9. Output matrices
# ============================================================

write.table(
  tibble(GeneID = rownames(rna_group_mean)) %>%
    bind_cols(as.data.frame(rna_group_mean, check.names = FALSE)),
  file = file.path(opt$outdir, "RNA_group_mean_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tibble(GeneID = rownames(ribo_group_mean)) %>%
    bind_cols(as.data.frame(ribo_group_mean, check.names = FALSE)),
  file = file.path(opt$outdir, "Ribo_group_mean_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tibble(GeneID = rownames(te_matrix)) %>%
    bind_cols(as.data.frame(te_matrix, check.names = FALSE)),
  file = file.path(opt$outdir, "TE_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tibble(GeneID = rownames(log2te_matrix)) %>%
    bind_cols(as.data.frame(log2te_matrix, check.names = FALSE)),
  file = file.path(opt$outdir, "log2TE_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 10. Summary table
# ============================================================

summary_table <- map_dfr(
  common_conditions,
  function(cond) {
    tibble(
      Condition = cond,
      total_common_genes = length(common_genes),
      genes_with_valid_TE = sum(!is.na(te_matrix[, cond])),
      genes_with_valid_log2TE = sum(!is.na(log2te_matrix[, cond]))
    )
  }
)

write.table(
  summary_table,
  file = file.path(opt$outdir, "TE_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 11. Prepare long-format data for plotting
# ============================================================

plot_dat <- log2te_matrix %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("GeneID") %>%
  pivot_longer(
    cols = -GeneID,
    names_to = "Condition",
    values_to = "log2TE"
  ) %>%
  filter(!is.na(log2TE)) %>%
  mutate(
    Condition = factor(Condition, levels = common_conditions)
  )

if (nrow(plot_dat) == 0) {
  warning("No valid log2TE values for plotting. Density plot and boxplot will not be generated.")
} else {
  
  # ============================================================
  # 12. Density plot
  # ============================================================
  
  p_density <- ggplot(
    plot_dat,
    aes(
      x = log2TE,
      color = Condition,
      #fill = Condition
    )
  ) +
    geom_density(
      #alpha = 0.20,
      linewidth = 0.8
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    ) +
    labs(
      title = "log2TE density plot",
      x = "log2(TE)",
      y = "Density",
      color = "Condition",
      fill = "Condition"
    )
  
  ggsave(
    filename = file.path(opt$outdir, "log2TE_density.pdf"),
    plot = p_density,
    width = 7,
    height = 5
  )
  
  
  # ============================================================
  # 13. Boxplot
  # ============================================================
  
  p_box <- ggplot(
    plot_dat,
    aes(
      x = Condition,
      y = log2TE,
      fill = Condition
    )
  ) +
    geom_boxplot(
      width = 0.65,
      outlier.size = 0.6,
      outlier.alpha = 0.5
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "none"
    ) +
    labs(
      title = "log2TE boxplot",
      x = "Condition",
      y = "log2(TE)"
    )
  
  ggsave(
    filename = file.path(opt$outdir, "log2TE_boxplot.pdf"),
    plot = p_box,
    width = 7,
    height = 7
  )
}

message("Done. Results were saved to: ", opt$outdir)