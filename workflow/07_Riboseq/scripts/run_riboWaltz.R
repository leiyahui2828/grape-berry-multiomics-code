#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(riboWaltz)
})

# -----------------------------
# 0. 辅助函数
# -----------------------------
get_ribowaltz_plot <- function(result, preferred_name = NULL) {
  candidate_names <- unique(c(preferred_name, "plot", "plot_Average"))
  for (nm in candidate_names) {
    if (!is.null(nm) && nm %in% names(result) && inherits(result[[nm]], "ggplot")) return(result[[nm]])
  }
  gg_names <- names(result)[vapply(result, function(x) inherits(x, "ggplot"), logical(1))]
  if (length(gg_names) > 0) return(result[[gg_names[1]]])
  stop("No ggplot object was found in the riboWaltz result.")
}

save_single_plot <- function(file, plot, width = 10, height = 7, plot_adjust = NULL) {
  if (!is.null(plot_adjust)) plot <- plot_adjust(plot)
  pdf(file = file, width = width, height = height)
  print(plot)
  dev.off()
}

save_sample_multipage_pdf <- function(file, sample_names, plot_function, width = 10, height = 7, plot_adjust = NULL) {
  pdf(file = file, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  for (s in sample_names) {
    result <- plot_function(s)
    p <- get_ribowaltz_plot(result, preferred_name = paste0("plot_", s)) + ggtitle(s)
    if (!is.null(plot_adjust)) p <- plot_adjust(p)
    print(p)
  }
}

write_ribowaltz_table <- function(result, file) {
  table_names <- c("count_dt", "dt", "plot_dt")
  for (nm in table_names) {
    if (nm %in% names(result)) {
      fwrite(as.data.table(result[[nm]]), file = file, sep = "\t")
      return(invisible(TRUE))
    }
  }
  invisible(FALSE)
}

# -----------------------------
# 1. 参数设置
# -----------------------------
option_list <- list(
  make_option(
    c("-b", "--bam_dir"),
    type = "character",
    help = "Input BAM directory. STAR Aligned.toTranscriptome.out.bam files are recommended."
  ),
  make_option(
    c("-g", "--gtf"),
    type = "character",
    help = "Reference annotation file in GTF or GFF format. It should match the annotation used for STAR index construction."
  ),
  make_option(
    c("-s", "--sample_info"),
    type = "character",
    help = "Sample information file without header. Column 1: group, column 2: sample."
  ),
  make_option(
    c("-f", "--fasta"),
    type = "character",
    default = NA,
    help = "Optional transcripts FASTA file. If provided, codon usage analysis will be performed using riboWaltz::codon_usage_psite()."
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = "riboWaltz_result",
    help = "Output directory [default: %default]"
  ),
  make_option(
    c("-t", "--threads"),
    type = "integer",
    default = 1,
    help = "Number of threads used by R, data.table, BLAS/OpenMP [default: %default]"
  ),
  make_option(
    c("--length_filter_mode"),
    type = "character",
    default = "custom",
    help = "Read length filtering mode: periodicity or custom [default: %default]"
  ),
  make_option(
    c("--periodicity_threshold"),
    type = "integer",
    default = 50,
    help = "Periodicity threshold used when length_filter_mode is periodicity [default: %default]"
  ),
  make_option(
    c("--min_len"),
    type = "integer",
    default = 25,
    help = "Minimum RPF length to keep [default: %default]"
  ),
  make_option(
    c("--max_len"),
    type = "integer",
    default = 35,
    help = "Maximum RPF length to keep [default: %default]"
  ),
  make_option(
    c("--flanking"),
    type = "integer",
    default = 6,
    help = "Minimum flanking length required around start/stop codons for psite() offset estimation [default: %default]"
  ),
  make_option(
    c("--extremity"),
    type = "character",
    default = "auto",
    help = "Read extremity used for P-site offset estimation: 5end, 3end, or auto [default: %default]"
  )

)

opt <- parse_args(OptionParser(option_list = option_list))

# -----------------------------
# 1.1 R 线程控制
# -----------------------------
if (is.null(opt$threads) || is.na(opt$threads) || opt$threads < 1) {
  stop("--threads must be a positive integer.")
}

THREADS <- as.integer(opt$threads)

Sys.setenv(
  OMP_NUM_THREADS = THREADS,
  OPENBLAS_NUM_THREADS = THREADS,
  MKL_NUM_THREADS = THREADS,
  VECLIB_MAXIMUM_THREADS = THREADS,
  NUMEXPR_NUM_THREADS = THREADS
)

if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::setDTthreads(THREADS)
}

# 检查必要参数
if (is.null(opt$bam_dir) || is.null(opt$gtf) || is.null(opt$sample_info)) {
  stop("--bam_dir, --gtf and --sample_info must be specified.")
}

# 创建输出目录
.dir_create <- function(x) dir.create(x, showWarnings = FALSE, recursive = TRUE)
.dir_create(opt$outdir)
.dir_create(file.path(opt$outdir, "tables"))
.dir_create(file.path(opt$outdir, "plots"))
#.dir_create(file.path(opt$outdir, "rdata"))

message(">>> Input BAM directory: ", opt$bam_dir)
message(">>> Input annotation file: ", opt$gtf)
message(">>> Input sample information file: ", opt$sample_info)
message(">>> Output directory: ", opt$outdir)

# -----------------------------
# 2. 检查 BAM 文件
# -----------------------------
bam_files <- list.files(opt$bam_dir, pattern = "\\.bam$", full.names = TRUE)

if (length(bam_files) == 0) {
  stop("No .bam files were found in the input BAM directory.")
}

message(">>> Number of BAM files detected: ", length(bam_files))
message(paste(basename(bam_files), collapse = "\n"))

# -----------------------------
# 3. 构建 riboWaltz 注释文件
# -----------------------------
message(">>> Step 1: Building riboWaltz annotation.")

annotation_dt <- create_annotation(gtfpath = opt$gtf)

fwrite(
  annotation_dt,
  file = file.path(opt$outdir, "tables", "riboWaltz_annotation.tsv"),
  sep = "\t"
)

#saveRDS(
#  annotation_dt,
#  file = file.path(opt$outdir, "rdata", "riboWaltz_annotation.rds")
#)

message(">>> Number of annotated transcripts: ", nrow(annotation_dt))

# -----------------------------
# 4. 读取 BAM 文件并转换为 reads list
# -----------------------------
message(">>> Step 2: Reading BAM files and converting them into a read list.")

reads_list <- bamtolist(
  bamfolder = opt$bam_dir,
  annotation = annotation_dt,
  transcript_align = TRUE
)

# 读取并检查样品分组信息；后续绘图直接使用 riboWaltz 内置 sample/multisamples 机制
sample_info <- fread(opt$sample_info, header = FALSE, data.table = TRUE)[, 1:2]
setnames(sample_info, c("group", "sample"))
sample_info[, `:=`(group = as.character(group), sample = as.character(sample))]

if (anyDuplicated(sample_info$sample)) {
  stop("Duplicated sample names in --sample_info: ",
       paste(unique(sample_info$sample[duplicated(sample_info$sample)]), collapse = ", "))
}
if (!setequal(sample_info$sample, names(reads_list))) {
  stop("Sample names in --sample_info do not match BAM-derived sample names.")
}

sample_names <- sample_info$sample
sample_groups <- split(sample_info$sample, sample_info$group)
reads_list <- reads_list[sample_names]

fwrite(
  sample_info,
  file = file.path(opt$outdir, "tables", "sample_info.used.tsv"),
  sep = "\t"
)

#saveRDS(
#  reads_list,
#  file = file.path(opt$outdir, "rdata", "reads_list.raw.rds")
#)

read_count_dt <- data.table(
  group = sample_info$group,
  sample = sample_names,
  reads = sapply(reads_list, nrow)
)

fwrite(
  read_count_dt,
  file = file.path(opt$outdir, "tables", "reads_count.raw.tsv"),
  sep = "\t"
)

message(">>> Raw read count summary has been generated.")

# -----------------------------
# 5. 统计 RPF 长度分布
# -----------------------------
message(">>> Step 3: Plotting RPF length distribution.")

# 单样品：多页 PDF，每页一个样品
save_sample_multipage_pdf(
  file = file.path(opt$outdir, "plots", "01_read_length_distribution.raw.sample.pdf"),
  sample_names = sample_names,
  width = 9,
  height = 7,
  plot_function = function(s) {
    rlength_distr(reads_list, sample = s, multisamples = "independent", plot_style = "dodge", cl = 99)
  }
)

# 分组：使用 riboWaltz 内置 list + multisamples="average" 处理重复
rlength_raw_group <- rlength_distr(
  reads_list,
  sample = sample_groups,
  multisamples = "average",
  plot_style = "dodge",
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "01_read_length_distribution.raw.group.pdf"),
  plot = get_ribowaltz_plot(rlength_raw_group) + theme_classic(),
  width = 9,
  height = 7
)

# 合并：所有样品作为一个整体平均
rlength_raw_merged <- rlength_distr(
  reads_list,
  sample = sample_names,
  multisamples = "average",
  plot_style = "dodge",
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "01_read_length_distribution.raw.merged.pdf"),
  plot = get_ribowaltz_plot(rlength_raw_merged) + theme_classic(),
  width = 9,
  height = 7
)

write_ribowaltz_table(
  rlength_distr(reads_list, sample = sample_names, multisamples = "independent", plot_style = "facet", cl = 99),
  file.path(opt$outdir, "tables", "read_length_distribution.raw.tsv")
)

# -----------------------------
# 6. 按 RPF 长度过滤 reads
# -----------------------------
message(">>> Step 4: Filtering RPFs by length: ", opt$length_filter_mode)

if (opt$length_filter_mode == "custom") {
  reads_list_filtered <- length_filter(
    data = reads_list,
    length_filter_mode = "custom",
    length_range = opt$min_len:opt$max_len
  )
} else if (opt$length_filter_mode == "periodicity") {
  reads_list_filtered <- length_filter(
    data = reads_list,
    length_filter_mode = "periodicity",
    periodicity_threshold = opt$periodicity_threshold
  )
} else {
  stop("Invalid --length_filter_mode. Please use 'custom' or 'periodicity'.")
}

#reads_list_filtered <- reads_list_filtered[sample_names]

#saveRDS(
#  reads_list_filtered,
#  file = file.path(opt$outdir, "rdata", "reads_list.filtered.rds")
#)

read_count_filtered_dt <- data.table(
  group = sample_info$group,
  sample = sample_names,
  reads = sapply(reads_list_filtered, nrow)
)

fwrite(
  read_count_filtered_dt,
  file = file.path(opt$outdir, "tables", "reads_count.filtered.tsv"),
  sep = "\t"
)

save_sample_multipage_pdf(
  file = file.path(opt$outdir, "plots", "02_read_length_distribution.filtered.sample.pdf"),
  sample_names = sample_names,
  width = 9,
  height = 7,
  plot_function = function(s) {
    rlength_distr(reads_list_filtered, sample = s, multisamples = "independent", plot_style = "dodge", cl = 99)
  }
)

rlength_filtered_group <- rlength_distr(
  reads_list_filtered,
  sample = sample_groups,
  multisamples = "average",
  plot_style = "dodge",
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "02_read_length_distribution.filtered.group.pdf"),
  plot = get_ribowaltz_plot(rlength_filtered_group) + theme_classic(),
  width = 9,
  height = 7
)

rlength_filtered_merged <- rlength_distr(
  reads_list_filtered,
  sample = sample_names,
  multisamples = "average",
  plot_style = "dodge",
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "02_read_length_distribution.filtered.merged.pdf"),
  plot = get_ribowaltz_plot(rlength_filtered_merged) + theme_classic(),
  width = 9,
  height = 7
)

write_ribowaltz_table(
  rlength_distr(reads_list_filtered, sample = sample_names, multisamples = "independent", plot_style = "facet", cl = 99),
  file.path(opt$outdir, "tables", "read_length_distribution.filtered.tsv")
)

# -----------------------------
# 7. 绘制 reads 5' 端和 3' 端分布热图
# -----------------------------
message(">>> Step 5: Plotting read 5' and 3' end distributions around start and stop codons.")

save_sample_multipage_pdf(
  file = file.path(opt$outdir, "plots", "03_read_extremity_heatmap.sample.pdf"),
  sample_names = sample_names,
  width = 12,
  height = 7,
  plot_function = function(s) {
    rends_heat(
      reads_list_filtered,
      annotation_dt,
      sample = s,
      multisamples = "independent",
      plot_style = "facet",
      cl = 99,
      utr5l = 25,
      cdsl = 40,
      utr3l = 25,
      colour = "blue"
    )
  }
)

ends_heat_group <- rends_heat(
  reads_list_filtered,
  annotation_dt,
  sample = sample_groups,
  multisamples = "average",
  plot_style = "facet",
  cl = 99,
  utr5l = 25,
  cdsl = 40,
  utr3l = 25,
  colour = "blue"
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "03_read_extremity_heatmap.group.pdf"),
  plot = get_ribowaltz_plot(ends_heat_group),
  width = 20,
  height = min(length(sample_groups) * 5, 20)
)

ends_heat_merged <- rends_heat(
  reads_list_filtered,
  annotation_dt,
  sample = sample_names,
  multisamples = "average",
  plot_style = "facet",
  cl = 99,
  utr5l = 25,
  cdsl = 40,
  utr3l = 25,
  colour = "blue"
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "03_read_extremity_heatmap.merged.pdf"),
  plot = get_ribowaltz_plot(ends_heat_merged),
  width = 12,
  height = 7
)

# -----------------------------
# 8. 计算 P-site offset
# -----------------------------
message(">>> Step 6: Estimating P-site offsets.")

psite_offset <- psite(
  reads_list_filtered,
  flanking = opt$flanking,
  extremity = opt$extremity,
  plot = TRUE,
  plot_dir = file.path(opt$outdir, "plots", "psite_offset_profiles")
)

fwrite(
  psite_offset,
  file = file.path(opt$outdir, "tables", "psite_offset.tsv"),
  sep = "\t"
)

#saveRDS(
#  psite_offset,
#  file = file.path(opt$outdir, "rdata", "psite_offset.rds")
#)

message(">>> P-site offset estimation has been completed.")

# -----------------------------
# 9. 为 reads 添加 P-site 位置信息
# -----------------------------
message(">>> Step 7: Adding P-site position information to reads.")

reads_psite_list <- psite_info(
  reads_list_filtered,
  psite_offset
)

#reads_psite_list <- reads_psite_list[sample_names]

#saveRDS(
#  reads_psite_list,
#  file = file.path(opt$outdir, "rdata", "reads_psite_list.rds")
#)

for (s in names(reads_psite_list)) {
  fwrite(
    as.data.table(reads_psite_list[[s]]),
    file = file.path(opt$outdir, "tables", paste0(s, ".reads_psite.tsv")),
    sep = "\t"
  )
}

# -----------------------------
# 10. 分析 P-site 在 5'UTR、CDS、3'UTR 中的分布
# -----------------------------
message(">>> Step 8: Analyzing P-site distribution across 5'UTR, CDS, and 3'UTR.")

save_sample_multipage_pdf(
  file = file.path(opt$outdir, "plots", "04_psite_region_distribution.sample.pdf"),
  sample_names = sample_names,
  width = 9,
  height = 7,
  plot_function = function(s) {
    region_psite(
      reads_psite_list,
      annotation_dt,
      sample = s,
      multisamples = "independent",
      plot_style = "stack",
      colour = c("lightblue", "lightgreen", "blue"),
      cl = 99
    )
  },
  plot_adjust = function(p) {
     p + theme(legend.text = element_text(margin = margin(r = 1)))
  }
)

# 堆叠图
region_result_group <- region_psite(
  reads_psite_list,
  annotation_dt,
  sample = sample_groups,
  multisamples = "average",
  plot_style = "stack", # dodge or stack
  colour = c("lightblue", "lightgreen", "blue"),
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "04_psite_region_distribution.stack.group.pdf"),
  plot = get_ribowaltz_plot(region_result_group) +
      theme(legend.text = element_text(margin = margin(r = 1))),
  width = min(length(sample_groups)*3 + 1, 20),
  height = 6
)

# 分组柱状图
region_result_group <- region_psite(
  reads_psite_list,
  annotation_dt,
  sample = sample_groups,
  multisamples = "average",
  plot_style = "dodge", # dodge or stack
  colour = c("lightblue", "lightgreen", "blue"),
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "04_psite_region_distribution.dodge.group.pdf"),
  plot = get_ribowaltz_plot(region_result_group) +
    theme(legend.text = element_text(margin = margin(r = 1))),
  width = min(length(sample_groups)*3 + 1, 20),
  height = 6
)

# 堆叠图
region_result_merged <- region_psite(
  reads_psite_list,
  annotation_dt,
  sample = sample_names,
  multisamples = "average",
  plot_style = "stack",
  colour = c("lightblue", "lightgreen", "blue"),
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "04_psite_region_distribution.stack.merged.pdf"),
  plot = get_ribowaltz_plot(region_result_merged) +
    theme(legend.text = element_text(margin = margin(r = 1))),
  width = 10,
  height = 6
)

# 分组柱状图
region_result_merged <- region_psite(
  reads_psite_list,
  annotation_dt,
  sample = sample_names,
  multisamples = "average",
  plot_style = "dodge",
  colour = c("lightblue", "lightgreen", "blue"),
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "04_psite_region_distribution.dodge.merged.pdf"),
  plot = get_ribowaltz_plot(region_result_merged) +
    theme(legend.text = element_text(margin = margin(r = 1))),
  width = 10,
  height = 6
)

write_ribowaltz_table(region_result_group, file.path(opt$outdir, "tables", "psite_region_distribution.tsv"))


# -----------------------------
# 11. 按 RPF 长度分析三核苷酸周期性
# -----------------------------
message(">>> Step 9: Analyzing three-nucleotide periodicity.")

save_sample_multipage_pdf(
  file = file.path(opt$outdir, "plots", "05_frame_distribution_by_length.sample.pdf"),
  sample_names = sample_names,
  width = 9,
  height = 7,
  plot_function = function(s) {
    frame_psite_length(
      reads_psite_list,
      annotation_dt,
      sample = s,
      multisamples = "independent",
      plot_style = "facet",
      region = "all",
      cl = 99
    )
  }
)

save_sample_multipage_pdf(
  file = file.path(opt$outdir, "plots", "05_frame_distribution_by_length_cds.sample.pdf"),
  sample_names = sample_names,
  width = 9,
  height = 7,
  plot_function = function(s) {
    frame_psite_length(
      reads_psite_list,
      annotation_dt,
      sample = s,
      multisamples = "independent",
      plot_style = "facet",
      region = "cds",
      cl = 99
    )
  }
)


frame_length_result_group <- frame_psite_length(
  reads_psite_list,
  annotation_dt,
  sample = sample_groups,
  multisamples = "average",
  plot_style = "facet",
  region = "all",
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "05_frame_distribution_by_length.group.pdf"),
  plot = get_ribowaltz_plot(frame_length_result_group),
  width = 9,
  height = length(sample_groups)*2 + 1
)

frame_length_result_merged <- frame_psite_length(
  reads_psite_list,
  annotation_dt,
  sample = sample_names,
  multisamples = "average",
  plot_style = "facet",
  region = "all",
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "05_frame_distribution_by_length.merged.pdf"),
  plot = get_ribowaltz_plot(frame_length_result_merged),
  width = 9,
  height = 6
)

write_ribowaltz_table(frame_length_result_merged, file.path(opt$outdir, "tables", "frame_distribution_by_length.tsv"))

# -----------------------------
# 12. 分析 CDS 区域整体 reading frame 分布
# -----------------------------
message(">>> Step 10: Plotting overall P-site frame distribution in CDS.")

save_sample_multipage_pdf(
  file = file.path(opt$outdir, "plots", "06_frame_distribution.cds.sample.pdf"),
  sample_names = sample_names,
  width = 8,
  height = 7,
  plot_function = function(s) {
    frame_psite(
      reads_psite_list,
      annotation_dt,
      sample = s,
      multisamples = "independent",
      plot_style = "facet",
      region = "cds",
      cl = 99
    )
  }
)

frame_result_group <- frame_psite(
  reads_psite_list,
  annotation_dt,
  sample = sample_groups,
  multisamples = "average",
  plot_style = "facet",
  region = "cds",
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "06_frame_distribution.cds.group.pdf"),
  plot = get_ribowaltz_plot(frame_result_group),
  width = 8,
  height = min(length(sample_groups) * 2 + 1, 20)
)

frame_result_merged <- frame_psite(
  reads_psite_list,
  annotation_dt,
  sample = sample_names,
  multisamples = "average",
  plot_style = "facet",
  region = "cds",
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "06_frame_distribution.cds.merged.pdf"),
  plot = get_ribowaltz_plot(frame_result_merged),
  width = 8,
  height = 6
)

write_ribowaltz_table(frame_result_merged, file.path(opt$outdir, "tables", "frame_distribution.cds.tsv"))

# -----------------------------
# 13. 绘制起始密码子和终止密码子附近的 metagene profile
# -----------------------------
message(">>> Step 11: Plotting P-site metagene profiles around start and stop codons.")

save_sample_multipage_pdf(
  file = file.path(opt$outdir, "plots", "07_psite_metaprofile_start_stop.sample.pdf"),
  sample_names = sample_names,
  width = 12,
  height = 7,
  plot_function = function(s) {
    metaprofile_psite(
      reads_psite_list,
      annotation_dt,
      sample = s,
      multisamples = "independent",
      plot_style = "facet",
      utr5l = 20,
      cdsl = 40,
      utr3l = 20,
      cl = 99
    )
  }
)

metaprofile_result_group <- metaprofile_psite(
  reads_psite_list,
  annotation_dt,
  sample = sample_groups,
  multisamples = "average",
  plot_style = "facet",
  utr5l = 20,
  cdsl = 40,
  utr3l = 20,
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "07_psite_metaprofile_start_stop.group.pdf"),
  plot = get_ribowaltz_plot(metaprofile_result_group),
  width = 20,
  height = min(length(sample_groups) * 5, 40)
)

metaprofile_result_group <- metaprofile_psite(
  reads_psite_list,
  annotation_dt,
  sample = sample_groups,
  multisamples = "average",
  plot_style = "overlap",
  utr5l = 20,
  cdsl = 40,
  utr3l = 20,
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "07_psite_metaprofile_start_stop.group_overlap.pdf"),
  plot = get_ribowaltz_plot(metaprofile_result_group),
  width = 20,
  height = 8
)

metaprofile_result_merged <- metaprofile_psite(
  reads_psite_list,
  annotation_dt,
  sample = sample_names,
  multisamples = "average",
  plot_style = "overlap",
  utr5l = 20,
  cdsl = 40,
  utr3l = 20,
  cl = 99
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "07_psite_metaprofile_start_stop.merged.pdf"),
  plot = get_ribowaltz_plot(metaprofile_result_merged),
  width = 20,
  height = 8
)


write_ribowaltz_table(metaprofile_result_merged, file.path(opt$outdir, "tables", "psite_metaprofile.tsv"))

save_sample_multipage_pdf(
  file = file.path(opt$outdir, "plots", "07_psite_metaprofile_start_stop.heatmap.sample.pdf"),
  sample_names = sample_names,
  width = 20,
  height = 7,
  plot_function = function(s) {
    metaheatmap_psite(
      reads_psite_list,
      annotation_dt,
      sample = s,
      multisamples = "independent",
      utr5l = 20,
      cdsl = 40,
      utr3l = 20,
      colour = "#333f50"
    )
  }
)

metaheatmap_result_group <- metaheatmap_psite(
  reads_psite_list,
  annotation_dt,
  sample = sample_groups,
  multisamples = "average",
  utr5l = 20,
  cdsl = 40,
  utr3l = 20,
  colour = "#333f50"
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "07_psite_metaprofile_start_stop.heatmap.group.pdf"),
  plot = get_ribowaltz_plot(metaheatmap_result_group),
  width = 20,
  height = min(length(sample_groups) * 2, 40)
)

metaheatmap_result_merged <- metaheatmap_psite(
  reads_psite_list,
  annotation_dt,
  sample = sample_names,
  multisamples = "average",
  utr5l = 20,
  cdsl = 40,
  utr3l = 20,
  colour = "#333f50"
)
save_single_plot(
  file = file.path(opt$outdir, "plots", "07_psite_metaprofile_start_stop.heatmap.merged.pdf"),
  plot = get_ribowaltz_plot(metaheatmap_result_merged),
  width = 20,
  height = 7
)


if (!is.na(opt$fasta)) {
  message(">>> Step 12: Analyzing codon usage based on P-site.")
  
  save_sample_multipage_pdf(
    file = file.path(opt$outdir, "plots", "08_codon_usage_psite.sample.pdf"),
    sample_names = sample_names,
    width = 20,
    height = 8,
    plot_function = function(s) {
      codon_usage_psite(
        data = reads_psite_list,
        annotation = annotation_dt,
        sample = s,
        multisamples = "independent",
        plot_style = "facet",
        site = "psite",
        frequency_normalization = TRUE,
        fastapath = opt$fasta,
        fasta_genome = FALSE,
        cl = 99,
        include_stop_codons = TRUE
      )
    },
    plot_adjust = function(p) {
      p + theme(legend.text = element_text(margin = margin(r = 1)))
    }
  )
  
  codon_usage_group <- codon_usage_psite(
    data = reads_psite_list,
    annotation = annotation_dt,
    sample = sample_groups,
    multisamples = "average",
    plot_style = "facet",
    site = "psite",
    frequency_normalization = TRUE,
    fastapath = opt$fasta,
    fasta_genome = FALSE,
    cl = 99,
    include_stop_codons = TRUE
  )
  
  #saveRDS(codon_usage_group, file = file.path(opt$outdir, "rdata", "codon_usage_psite.group.rds"))
  
  save_single_plot(
    file = file.path(opt$outdir, "plots", "08_codon_usage_psite.group.pdf"),
    plot = get_ribowaltz_plot(codon_usage_group) + theme(legend.text = element_text(margin = margin(r = 1))),
    width = 20,
    height = min(length(sample_groups) * 3, 40)
  )
 
  codon_usage_merged <- codon_usage_psite(
    data = reads_psite_list,
    annotation = annotation_dt,
    sample = sample_names,
    multisamples = "average",
    plot_style = "facet",
    site = "psite",
    frequency_normalization = TRUE,
    fastapath = opt$fasta,
    fasta_genome = FALSE,
    cl = 99,
    include_stop_codons = TRUE
  )

  save_single_plot(
    file = file.path(opt$outdir, "plots", "08_codon_usage_psite.merged.pdf"),
    plot = get_ribowaltz_plot(codon_usage_merged) + theme(legend.text = element_text(margin = margin(r = 1))),
    width = 20,
    height = 8
  )
  write_ribowaltz_table(codon_usage_merged, file.path(opt$outdir, "tables", "codon_usage_psite.tsv"))
  
  # pairwise group contrasts: save all codon usage scatter plots into one multi-page PDF
  if (length(sample_groups) >= 2) {
    message(">>> Plotting pairwise group codon usage scatter plots into one PDF.")
    contrast_pairs <- combn(names(sample_groups), 2, simplify = FALSE)
    pdf(file = file.path(opt$outdir, "plots", "08_codon_usage_psite.contrast.group.pdf"), width = 8, height = 7)
    for (contrast_pair in contrast_pairs) {
      pair_name <- paste(contrast_pair, collapse = "_vs_")
      pair_file_name <- gsub("[^A-Za-z0-9_.-]+", "_", pair_name)
      codon_usage_contrast <- codon_usage_psite(
        data = reads_psite_list,
        annotation = annotation_dt,
        sample = sample_groups,
        multisamples = "average",
        contrast_sample = contrast_pair,
        site = "psite",
        frequency_normalization = TRUE,
        fastapath = opt$fasta,
        fasta_genome = FALSE,
        cl = 99,
        include_stop_codons = TRUE,
        label_scatter = TRUE,
        label_number = 5
      )
 
      print(get_ribowaltz_plot(codon_usage_contrast) + theme_classic() + ggtitle(pair_name))
    }
    dev.off()
  } else {
    message(">>> Less than two groups were detected. Skipping pairwise codon usage contrast scatter plots.")
  }
  
  
  message(">>> Codon usage analysis has been completed.")
} else {
  message(">>> Step 12: No --fasta was provided. Skipping codon usage analysis.")
}

# -----------------------------
# 15. 输出 R 运行环境信息
# -----------------------------
sink(file.path(opt$outdir, "sessionInfo.txt"))
sessionInfo()
sink()

message(">>> riboWaltz analysis has been completed successfully.")
message(">>> Output directory: ", normalizePath(opt$outdir))
