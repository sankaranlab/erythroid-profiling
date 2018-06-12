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
mapvec <- make.unique(as.character(mapdf[,2])); names(mapvec) <- as.character(mapdf[,1])
colnames(RNA.counts) <- mapvec[as.character(colnames(RNA.counts))]

RNA.counts$genes <- make.unique(read.table("../../data/genes.tsv", stringsAsFactors = FALSE)[,2])

# Polish
qcRNA$PercentReadsInGenes <- round(as.numeric(as.character(qcRNA$ReadsAssignedToGenes))/as.numeric(as.character(qcRNA$TotalReads)) * 100, 1)
write.table(qcRNA, file = "../../downloads/RNAseq_QC.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

write.table(RNA.counts, file = "../../processed/RNAseq_rawGeneCounts.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
