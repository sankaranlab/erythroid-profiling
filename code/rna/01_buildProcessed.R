library(dplyr)
library(data.table)

# Handle QC first
lapply(list.files("../../data/RNA_logs/"), function(name){
  sample <- gsub("_RNALog.final.out", "", name)
  file <- paste0("../../data/RNA_logs/", name)
  dat.in <- readLines(file)
  dat.in <- dat.in[grep("[|]", dat.in, invert=FALSE)]
  t <- read.table(textConnection(paste(dat.in, collapse="\n")), header=FALSE, sep = "\t")
  
  data.frame( Sample = sample, TotalReads = t[5,2], UniquelyMapping = t[7,2])
  
}) %>% rbindlist() %>% data.frame(stringsAsFactors = FALSE) -> qcRNA

# Now deal with expression
rnafiles <- list.files("../../data/RNA_genes", pattern=".tab", full.names = TRUE)
samplenames <- gsub("_RNAReadsPerGene.out.tab", "", list.files("../../data/RNA_genes", pattern=".tab"))
RNA.counts <- sapply(rnafiles, function(file) {
  data.frame(fread(paste0("cat < ", file), skip = 4))[,2]
}) %>% data.matrix()
colnames(RNA.counts) <- samplenames
RNA.counts <- data.frame(RNA.counts)
qcRNA$ReadsAssignedToGenes <- colSums(RNA.counts)[qcRNA$Sample]

RNA.counts$genes <- make.unique(read.table("../../data/genes.tsv", stringsAsFactors = FALSE)[,2])

# Polish
qcRNA$PercentReadsInGenes <- round(as.numeric(as.character(qcRNA$ReadsAssignedToGenes))/as.numeric(as.character(qcRNA$TotalReads)) * 100, 1)
write.table(qcRNA, file = "../../downloads/RNAseq_QC.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

write.table(RNA.counts, file = "../../processed/RNAseq_rawGeneCounts.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
