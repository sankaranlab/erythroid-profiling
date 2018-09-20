library(dplyr)
library(BuenColors)
library(data.table)

# Import raw expression values
raw <- read.table("../../processed/RNAseq_rawGeneCounts.tsv", header = TRUE)
RNA.counts <- raw[,1:28]
meta <- stringr::str_split_fixed(colnames(RNA.counts), "_", 4)
genes <- as.character(raw[,29])
cpm <- sweep(RNA.counts, 2, colSums(RNA.counts), FUN="/") * 1000000
log2cpm <- log2(cpm + 1)

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")


# Population mean and SE
pops <- paste0("P", as.character(1:8))
lapply(pops, function(i){
  n <- sum(meta[,3] == i)
  means <- rowMeans(log2cpm[,meta[,3] == i])
  ses <- sqrt(matrixStats::rowVars(data.matrix(log2cpm[,meta[,3] == i])))/sqrt(n)
  data.frame(Population = i, 
             Gene = genes,
             Mean = means, 
             SEM = ses)
}) %>% rbindlist() %>% data.frame() -> odf

# Import DE results
lapply(list.files("../../downloads/RNA_DESeq2/", full.names = TRUE), function(file){
  dt <- fread(file) %>% data.frame() %>% mutate(stat = -1*log10(padj) * log2FoldChange) %>% arrange(stat)
  c(head(dt$gene, 10), tail(dt$gene, 10))
}) %>% unlist() %>% unique() -> signature

odf[,c("Population", "Gene", "Mean")] %>% filter(Gene %in% signature) %>% tidyr::spread(key = Population, value = Mean) -> signatureMatrix
write.table(signatureMatrix, file = "../../processed/erythroid_signature_matrix.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

