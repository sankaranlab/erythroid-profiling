library(dplyr)
library(data.table)

files <- list.files("../../data/ATAC_data/sampleQC/")

makeTableATAC <- function(file){
  t <- read.table(paste0("../../data/ATAC_data/sampleQC/", file), header = TRUE)
  ss <- stringr::str_split_fixed(t$Sample, "_", 5)[1,]
  t$Sample <- paste0(ss[3], "-", ss[4])
  t
}

lapply(files, makeTableATAC) %>% rbindlist() %>% data.frame() %>% arrange(Sample) -> df
write.table(df, file = "../../processed/ATAC-QC-TableS2.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
