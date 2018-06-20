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

cpm_ct <- sapply(unique(unname(celltype)), function(cell){
  x <- rowSums(data.matrix(RNA.counts.all)[,cell == celltype], na.rm = TRUE)
  x/sum(x) * 1000000
})
cpm_ct <- data.frame(round(cpm_ct, 2))
cpm_ct$genes <- RNA.counts.all$genes

write.table(cpm_ct[,c("HSC", "MPP", "CMP", "MEP", paste0("P", seq(1:8)), "genes")], 
            file = "../../processed/expression_CPM_erythroid.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
