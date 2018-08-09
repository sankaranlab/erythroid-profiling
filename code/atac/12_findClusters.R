library(BuenColors)
library(dplyr)
library(DESeq2)
library(BiocParallel)
library(data.table)
library(annotables)
library(circlize)
library(ComplexHeatmap)

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

scaleRows <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
}

ATAC.counts.filt.Z <- scaleRows(ATAC.counts.df[boo,])
peaks.filt <- peaks[boo,]

if(FALSE){
  findBest <- cluster::clusGap(ATAC.counts.filt.Z, FUN = kmeans, K.max = 10, B = 3, nstart = 2)
  qplot(1:10,findBest$Tab[,3])
  
  df <- data.frame(K = 1:10, Gap = findBest$Tab[,3])
  #saveRDS(df, file = "../../processed/kmean-scree-atac.rds")
  df <- readRDS("../../processed/kmean-scree-atac.rds")
  p1 <- ggplot(df, aes(x = K, y = Gap)) + geom_point() +
    pretty_plot(fontsize = 8) + L_border() + labs(x = "K-means", y = "Gap Statistic") +
    geom_vline(xintercept = 7, linetype = 2) +
    scale_x_continuous(breaks = c(2,4,6,8,10,12,14))
  
  cowplot::ggsave(p1, file = "../../plots/atacseq_kmeans_scree.pdf", width = 2.2, height =2)
}


km <- kmeans(ATAC.counts.filt.Z, centers = 7, nstart = 2)
vals <- c("K5","K7","K6", "K2","K3","K4","K1")
km.cluster <- factor(vals[as.numeric(km$cluster)], levels = sort(vals))

# Sort DF
pops_order <- c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8)))
ATAC.condition <- factor(pops_order, levels = pops_order)
sdf <- data.frame(rep = pops_order,
                  pop = ATAC.condition) %>% arrange(pop)

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")
all_color_maps <- c(ejc_color_maps[c("HSC", "MPP", "CMP")], "MEP" = "#FF81AF" ,eryth_color_maps)

ha_col <- HeatmapAnnotation(cell = pull(sdf, pop), col = list(cell = all_color_maps))
ATAC.counts.filt.Z.plot <- ATAC.counts.filt.Z
ATAC.counts.filt.Z.plot[ATAC.counts.filt.Z.plot > 4] <- 4
ATAC.counts.filt.Z.plot[ATAC.counts.filt.Z.plot < -4] <- -4


pdf(file="../../plots/atacseq-clustering_allpops.pdf", width = 5, height = 7)  
par(cex.main=0.8,mar=c(1,1,1,1))
hm <- Heatmap(ATAC.counts.filt.Z.plot[,pull(sdf, rep)],
              col=as.character(jdb_palette("solar_extra",type="continuous")),
        cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 0),
        column_names_gp = gpar(fontsize = 0),
        top_annotation = ha_col,
        split = km.cluster, show_heatmap_legend = FALSE,
        name = "")
hm
dev.off()

if(FALSE){
  write.table(data.frame(data.frame(peaks.filt)[,c(1,2,3)], 
                         cluster = km.cluster),
              file = "../../processed/ATACseq-Kmeans7-clusterID-allpops.tsv", 
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}


