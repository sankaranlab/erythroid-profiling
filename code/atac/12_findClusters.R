library(BuenColors)
library(dplyr)
library(DESeq2)
library(BiocParallel)
library(data.table)
library(annotables)

register(MulticoreParam(2))

# Import raw expression values
ATAC.counts <- read.table("../../data/ATAC_data/combined_12.counts.tsv", header = TRUE)
peaks <- data.frame(fread("../../data/ATAC_data/ery_only.bed"))

# RNA DEseq2 setup
ATAC.counts.df <- as.data.frame(log2(sweep(data.matrix(ATAC.counts), 2, colSums(data.matrix(ATAC.counts)), FUN="/") * 1000000 + 1))

# Establish column data
ATAC.condition <- colnames(ATAC.counts)
colData <- as.data.frame(ATAC.condition)
row.names(colData) <- colnames(ATAC.counts.df)

# Import differentially accessible peaks
lapply(list.files("../../downloads/ATAC_DESeq2/", full.names = TRUE), function(x){
  data.frame(fread(x, header = TRUE, stringsAsFactors = FALSE))[,c(1,2,3)]
}) %>% rbindlist() %>% data.frame() -> diffPeakDf

boo <- peaks$V2 %in% diffPeakDf$start

ATAC.counts.filt <- ATAC.counts.df[boo,]
peaks.filt <- peaks[boo,]

if(TRUE){
  findBest <- cluster::clusGap(ATAC.counts.filt, FUN = kmeans, K.max = 15, B = 3, nstart = 2)
  qplot(1:10,findBest$Tab[,3])
  
  df <- data.frame(K = 1:10, Gap = findBest$Tab[,3])
  saveRS(df, file = "../../processed/kmean-scree-atac.rds")
  p1 <- ggplot(df, aes(x = K, y = Gap)) + geom_point() +
    pretty_plot(fontsize = 8) + L_border() + labs(x = "K-means", y = "Gap Statistic") +
    geom_vline(xintercept = 7, linetype = 2) +
    scale_x_continuous(breaks = c(2,4,6,8,10))
  
  cowplot::ggsave(p1, file = "../../plots/atacseq_kmeans_scree.pdf", width = 2.2, height =2)
}


