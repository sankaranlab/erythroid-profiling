library(BuenColors)
library(dplyr)
library(data.table)
library(annotables)

# Filter sex chromosomes
mdf <- merge(read.table("../../data/genes.tsv", header = FALSE), grch37, by.x = "V1", by.y = "ensgene")
symbol_keep <- unique(mdf[mdf$chr %in% c(as.character(1:22), "X"),"V2"])

# Import raw expression values for HSC-P8
load("../../processed/serialized/RNAseq_wCorces.rda")
RNA.counts.all <- RNA.counts.all[RNA.counts.all$genes %in% symbol_keep,]
genes <- RNA.counts.all$genes

# Collapse multiple samples into each cell type
test <- melt(RNA.counts.all) %>% as.data.table
test$ID <- celltype[test$variable]
test %>%
  group_by(genes, ID) %>% 
  summarise(allcounts = sum(value)) -> out
RNA.counts <- dcast(out,genes~ID)

if (FALSE){
  # Import raw expression values
  raw <- read.table("../../processed/RNAseq_rawGeneCounts.tsv", header = TRUE)
  raw <- raw[raw[,29] %in% symbol_keep,]
  RNA.counts <- raw[,1:28]
  meta <- stringr::str_split_fixed(colnames(RNA.counts), "_", 4)
  genes <- raw[,29]
  
  # Establish column data
  RNA.condition <- meta[,3]
  colData <- as.data.frame(RNA.condition)
  row.names(colData) <- colnames(RNA.counts.df)
}

# RNA DEseq2 setup
genes <- RNA.counts$genes
RNA.counts.df <- as.data.frame(log2(sweep(data.matrix(RNA.counts[,-1]), 2, colSums(data.matrix(RNA.counts[,-1])), FUN="/") * 1000000 + 1))


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
    geom_vline(xintercept = 7, linetype = 2) +
    scale_x_continuous(breaks = c(2,4,6,8,10))
  
  cowplot::ggsave(p1, file = "../../plots/rnaseq_kmeans_scree_allpops.pdf", width = 2.2, height =2)
}

km <- kmeans(RNA.counts.filt.Z, centers = 7, nstart = 100)
vals <- c("K7","K6","K2", "K4","K3","K1","K5")
km.cluster <- factor(vals[as.numeric(km$cluster)], levels = sort(vals))

# Sort DF
pops_order <- c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8)))
RNA.condition <- factor(pops_order, levels = pops_order)
sdf <- data.frame(rep = pops_order,
                  pop = RNA.condition) %>% arrange(pop)

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")
all_color_maps <- c(ejc_color_maps[c("HSC", "MPP", "CMP")], "MEP" = "#FF81AF" ,eryth_color_maps)

ha_col <- HeatmapAnnotation(cell = pull(sdf, pop), col = list(cell = all_color_maps))
RNA.counts.filt.Z.plot <- RNA.counts.filt.Z
RNA.counts.filt.Z.plot[RNA.counts.filt.Z.plot > 4] <- 4
RNA.counts.filt.Z.plot[RNA.counts.filt.Z.plot < -4] <- -4


pdf(file="../../plots/rnaseq-clustering_allpops.pdf", width = 5, height = 7)  
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
              file = "../../processed/RNAseq-Kmeans7-clusterID-allpops.tsv", 
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}
pops_order <- c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8)))
head(RNA.counts.filt[,pops_order])

all_df <- data.frame( read.table("../../processed/RNAseq-Kmeans7-clusterID-allpops.tsv", header = TRUE),
                      round(RNA.counts.filt[,pops_order], 2))

write.table(all_df,
            file = "../../processed/RNA-seq_suppTable.tsv", 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
