library(dplyr)
library(data.table)

# Now deal with expression
rnafiles <- list.files("data", pattern=".tab", full.names = TRUE)
samplenames <- gsub("_RNAReadsPerGene.out.tab", "", list.files("data", pattern=".tab"))
RNA.counts <- sapply(rnafiles, function(file) {
  data.frame(fread(paste0("cat < ", file), skip = 4))[,2]
}) %>% data.matrix()
colnames(RNA.counts) <- samplenames
RNA.counts <- data.frame(RNA.counts)
RNA.counts$genes <- make.unique(read.table("../../data/genes.tsv", stringsAsFactors = FALSE)[,2])
colnames(RNA.counts) <- gsub("ReadsPerGene.out.tab", "", colnames(RNA.counts))

write.table(RNA.counts, file = "output/TMCC2-RNAseq_rawGeneCounts.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
