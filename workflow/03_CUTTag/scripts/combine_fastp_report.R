suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(argparse)
})

# 创建命令行参数解析器
parser <- ArgumentParser()
parser$add_argument("--input", type="character", nargs='+', help="input fastp JSON reports")
parser$add_argument("--output", type="character", default="report", help="output folder to store fastp summary TSV and RData files")

# 解析命令行参数
args <- parser$parse_args()

# 获取输入文件和输出目录
json_files <- args$input
output_folder <- args$output

# 确保指定了输入文件
if (length(json_files) == 0) {
  stop("No input files provided.")
}

# 创建输出文件夹
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# 输出文件路径
output_file <- file.path(output_folder, "fastp_summary.tsv")

# 将要提取的字段列在这里
fields <- c(
  "total_reads",
  "total_bases",
  "q30_rate",
  "gc_content"
)

# 初始化一个空 data.frame 来存储结果
fastp_summary <- data.frame()

# 遍历所有 JSON 文件
for (json_file in json_files) {
  # 读取 JSON 文件
  data <- fromJSON(json_file)
  
  # 从文件名中提取样本名称
  sample_name <- gsub("_fastp.json", "", basename(json_file))
  
  # 提取过滤前和过滤后的所需字段
  before_filtering <- data$summary$before_filtering[fields]
  after_filtering <- data$summary$after_filtering[fields]
  
  # 将结果添加到 data.frame
  result <- data.frame(sample = sample_name)
  result <- cbind(result, before_filtering, after_filtering)
  fastp_summary <- rbind(fastp_summary, result)
}

# 重命名列名
colnames(fastp_summary) <- c("sample",
                             paste0("before_", fields),
                             paste0("after_", fields))

# 排序
fastp_summary <- fastp_summary %>% arrange(sample)

# 将结果写入指定的 TSV 文件
write.table(fastp_summary, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

# 存储 fastp_summary 数据到 RData 文件中
save(fastp_summary, file = file.path(output_folder, "fastp_summary.RData"))
