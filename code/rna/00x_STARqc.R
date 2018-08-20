library(data.table)
library(dplyr)

samples <- gsub("Log.final.out", "", list.files("../../data/RNA_logs/", pattern = "final.out$"))

lapply(samples, function(sample){
  log <- paste0("../../data/RNA_logs//", sample, "Log.final.out")
  lines <- readLines(log)
  avg_mapped_length <- strsplit(lines[[11]], "\t")[[1]][2]
  reads_in  <- strsplit(lines[[6]], "\t")[[1]][2]
  align_percent <-  as.numeric(gsub("%", "", strsplit(lines[[10]], "\t")[[1]][2])) * 0.01
  
  data.frame(sample, reads_in, avg_mapped_length, align_percent, aligned_reads = align_percent*as.numeric(reads_in))
  
}) %>% rbindlist() %>% data.frame() -> sumStats


