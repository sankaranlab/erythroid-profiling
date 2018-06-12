library(dplyr)
library(data.table)

# Now deal with expression
rnafiles <- list.files("../../data/corces/RNAseq_counts/", pattern=".tab", full.names = TRUE)
samplenames <- gsub("ReadsPerGene.out.tab", "", list.files("../../data/corces/RNAseq_counts", pattern=".tab"))
RNA.counts <- sapply(rnafiles, function(file) {
  data.frame(fread(paste0("cat < ", file), skip = 4))[,2]
}) %>% data.matrix()
colnames(RNA.counts) <- samplenames
RNA.counts <- data.frame(RNA.counts)

# Import mapping / translate SRR files
mapdf <- read.table("../../data/corces/RNAseq_Corces.txt", sep = "\t")
mapvec <- as.character(mapdf[,2]); names(mapvec) <- as.character(mapdf[,1])

df <- data.frame(fread("../../processed/RNAseq_rawGeneCounts.tsv"))

RNA.counts.all <- cbind(RNA.counts, df)
celltype <- c(mapvec[colnames(RNA.counts)], stringr::str_split_fixed(colnames(df), "_", 4)[1:28,3])

save(RNA.counts.all, celltype, file = "../../processed/serialized/RNAseq_wCorces.rda")
