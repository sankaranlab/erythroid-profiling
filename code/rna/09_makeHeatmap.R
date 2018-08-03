library(BuenColors)
library(dplyr)
library(data.table)
library(annotables)

# Filter sex chromosomes
mdf <- merge(read.table("../../data/genes.tsv", header = FALSE), grch37, by.x = "V1", by.y = "ensgene")
symbol_keep <- unique(mdf[mdf$chr %in% c(as.character(1:22), "X"),"V2"])

# Import raw expression values
raw <- read.table("../../processed/RNAseq_rawGeneCounts.tsv", header = TRUE)
raw <- raw[raw[,29] %in% symbol_keep,]
RNA.counts <- raw[,1:28]
meta <- stringr::str_split_fixed(colnames(RNA.counts), "_", 4)
genes <- raw[,29]

# RNA DEseq2 setup
RNA.counts.df <- as.data.frame(log2(sweep(data.matrix(RNA.counts), 2, colSums(data.matrix(RNA.counts)), FUN="/") * 1000000 + 1))

# Establish column data
RNA.condition <- meta[,3]
colData <- as.data.frame(RNA.condition)
row.names(colData) <- colnames(RNA.counts.df)

# Import differentially expressed genes
lapply(list.files("../../downloads/RNA_DESeq2", full.names = TRUE), function(x){
  data.frame(g = read.table(x, header = TRUE, stringsAsFactors = FALSE)[,1])
}) %>% rbindlist() %>% data.frame() -> diffGeneDf
differentialGenes <- unique(as.character(diffGeneDf[,1]))

boo <- genes %in% differentialGenes

RNA.counts.filt <- RNA.counts.df[boo,]
genes.filt <- genes[boo]

scaleRows <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
}

RNA.counts.filt.Z <- scaleRows(RNA.counts.filt)

if(FALSE){
  findBest <- cluster::clusGap(RNA.counts.filt.Z, FUN = kmeans, K.max = 10, B = 3, nstart = 2)
  qplot(1:10,findBest$Tab[,3])
  
  df <- data.frame(K = 1:10, Gap = findBest$Tab[,3])
  
  p1 <- ggplot(df, aes(x = K, y = Gap)) + geom_point() +
    pretty_plot(fontsize = 8) + L_border() + labs(x = "K-means", y = "Gap Statistic") +
    geom_vline(xintercept = 5, linetype = 2) +
    scale_x_continuous(breaks = c(2,4,6,8,10))
  
  cowplot::ggsave(p1, file = "../../plots/rnaseq_kmeans_scree.pdf", width = 2.2, height =2)
}

km <- kmeans(RNA.counts.filt.Z, centers = 5, nstart = 5)
vals <- c("K4", "K2", "K1", "K5", "K3")
km.cluster <- factor(vals[as.numeric(km$cluster)], levels = sort(vals))

# Sort DF
sdf <- data.frame(rep = colnames(RNA.counts.filt.Z),
                  pop = RNA.condition) %>% arrange(pop)

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")
ha_col <- HeatmapAnnotation(cell = pull(sdf, pop), col = list(cell = eryth_color_maps))
RNA.counts.filt.Z.plot <- RNA.counts.filt.Z
RNA.counts.filt.Z.plot[RNA.counts.filt.Z.plot > 4] <- 4
RNA.counts.filt.Z.plot[RNA.counts.filt.Z.plot < -4] <- -4


pdf(file="../../plots/rnaseq-clustering.pdf", width = 5, height = 7)  
par(cex.main=0.8,mar=c(1,1,1,1))
hm <- Heatmap(RNA.counts.filt.Z.plot[,pull(sdf, rep)],
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
  write.table(data.frame(gene = genes.filt, 
                         cluster = km.cluster),
              file = "../../processed/RNAseq-Kmeans5-clusterID.tsv", 
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}
