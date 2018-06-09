library(reshape2)
library(dplyr)


# Import RNA-seq
lapply(list.files("../../downloads/RNA_DESeq2"), function(x){
  tab <- read.table(paste0("../../downloads/RNA_DESeq2", "/", x), header = TRUE)
  
})


