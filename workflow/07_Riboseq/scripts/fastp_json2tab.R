#!/usr/bin/env Rscript

library(argparser, quietly=TRUE)# Create a parser
p <- arg_parser("Generate fastq stat infomation from fastp json files.")

# Add command line arguments
p <- add_argument(p, "json", help="input: json.list", type="character" )
p <- add_argument(p, "output", help="output file", type="character")

argv <- parse_args(p)

input <- argv$json
output <- argv$output

#input <- "./fastp.json.list"
#output <- "fastp.stat"

library("rjson")

json <- read.delim(input,header = F)

outtab <- data.frame()
for ( i in rownames(json)) { 
  name <- json[i,1]
  jfile <- json[i,2]
  res <-as.data.frame(fromJSON(file = jfile)$summary)
  row.names(res) <- name

  res <- res[, sapply(res, is.numeric), drop = FALSE]
  # 对数值列保留 4 位小数
  num_cols <- sapply(res, is.numeric)
  res[num_cols] <- lapply(res[num_cols], round, digits = 4)
  res <- cbind(Sample = name, res)
  outtab <- rbind(outtab, res)
  
}

options(scipen = 100)
write.table(outtab , file = output, quote = F, sep = "\t", row.names = F, col.names = T,)
