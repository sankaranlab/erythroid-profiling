library(dplyr)
library(GenomicRanges)
library(BiocParallel)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(BuenColors)
set.seed(14651)

# Import
bed <- data.frame(fread("../../data/corces/panHeme.bed", header = FALSE))
colnames(bed) <- c("chr", "start", "end")
peak_gr <- makeGRangesFromDataFrame(bed)

# Import raw expression values
raw <- data.frame(fread("../../data/corces/atac_combinedPopulations/panHeme.counts.tsv", header = TRUE))
ATAC.counts <- data.matrix(raw)
meta <- colnames(ATAC.counts)

# Normalize
cpm <- sweep(raw, 2, colSums(raw), FUN="/") * 1000000

# Import GATA1 peaks
gata1_gr <- diffloop::bedToGRanges("../../data/public-chipseq/GSM1067274_Erythroid_GATA1_peaks.bed")
ov1 <- findOverlaps(gata1_gr, peak_gr)

atac <- log2(cpm[1:dim(cpm)[1] %in% subjectHits(ov1), c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(2:8)))] +1)
atac.norm <- t(apply(atac, 1, function(x)(x-min(x))/(max(x)-min(x))))

# Figure out ideal numbers of k-means
findBest <- cluster::clusGap(atac.norm, FUN = kmeans, K.max = 10, B = 3, nstart = 100)
qplot(1:10,findBest$Tab[,3])

if(FALSE){
  df <- data.frame(K = 1:10, Gap = findBest$Tab[,3])
  
  p1 <- ggplot(df, aes(x = K, y = Gap)) + geom_point() +
    pretty_plot(fontsize = 8) + L_border() + labs(x = "K-means", y = "Gap Statistic") +
    geom_vline(xintercept = 9, linetype = 2)
  
  cowplot::ggsave(p1, file = "../../plots/gata1_peaks_kmeans_scree.pdf", width = 2.2, height =2)
}

km <- kmeans(atac.norm, centers = 4, nstart = 100)

groups <- c("g1", "g2", "g3", "g4")
all_clusters <- groups[as.numeric(km$cluster)]

big_df <- data.frame(bed[1:dim(cpm)[1] %in% subjectHits(ov1),], all_clusters, atac.norm)
refined_gr <- peak_gr[1:dim(cpm)[1] %in% subjectHits(ov1)]

# Import other co-factors
klf1_gr <- diffloop::bedToGRanges("../../data/public-chipseq/GSM1067275_Erythroid_KLF1_peaks.bed")
tal1_gr <- diffloop::bedToGRanges("../../data/public-chipseq/GSM1067277_Erythroid_TAL1_peaks.bed")
nfe2_gr <- diffloop::bedToGRanges("../../data/public-chipseq/GSM1067276_Erythroid_NFE2_peaks.bed")

big_df$nfe2  <- 1:dim(big_df)[1] %in% queryHits(findOverlaps(refined_gr, nfe2_gr))
big_df$tal1  <- 1:dim(big_df)[1] %in% queryHits(findOverlaps(refined_gr, tal1_gr))
big_df$klf1  <- 1:dim(big_df)[1] %in% queryHits(findOverlaps(refined_gr, klf1_gr))


big_df %>% group_by(all_clusters) %>%
  summarize(mHSC = mean(HSC),
            mMPP = mean(MPP),
            mCMP = mean(CMP), 
            mMEP = mean(MEP),
            mP2 = mean(P2),
            mP3 = mean(P3), 
            mP4 = mean(P4),
            mP5 = mean(P5), 
            mP6 = mean(P6), 
            mP7 = mean(P7), 
            mP8 = mean(P8)) -> summary_df

big_df %>% group_by(all_clusters) %>%
  summarise(pNFE2 = mean(nfe2), 
            pTAL1 = mean(tal1),
            pKLF1 = mean(klf1))

mm <- reshape2::melt(summary_df, id.vars = "all_clusters")

ggplot(mm, aes(x = variable, y = value, group = all_clusters, color = all_clusters)) +
  geom_point() + geom_line()
