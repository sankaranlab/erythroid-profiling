library(data.table)
library(dplyr)
library(stringr)

samples <- gsub("Log.final.out", "", list.files("../../data/RNA_logs/", pattern = "final.out$"))

lapply(samples, function(sample){
  log <- paste0("../../data/RNA_logs/", sample, "Log.final.out")
  lines <- readLines(log)
  avg_mapped_length <- strsplit(lines[[11]], "\t")[[1]][2]
  reads_in  <- strsplit(lines[[6]], "\t")[[1]][2]
  align_percent <-  as.numeric(gsub("%", "", strsplit(lines[[10]], "\t")[[1]][2])) * 0.01
  val <- paste0("../../data/RNA_genes/", sample, "ReadsPerGene.out.tab")
  quants <- data.frame(fread(val))
  n_noFeature <- (quants$V2[4])/sum(quants$V2[5:dim(quants)[1]])
  ss <- str_split_fixed(sample, "_", 6)
  sample2 <- paste0(ss[1,3], "-", ss[1,4])
  data.frame(sample2, reads_in, align_percent, percent_in_gene = 1 - n_noFeature)
  
}) %>% rbindlist() %>% data.frame() -> sumStats
