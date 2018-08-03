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

# Find overlaps with fine-mapped SNPS
import_to_df <- function(trait){
  file <- paste0("../../data/FMsnps/",trait,"_PP001_betas.bed")
  tab <- read.table(file)
  tab <- tab[tab$V5 > 0.1,]
  tab$trait <- trait
  colnames(tab) <- c("chr", "start", "end", "region", "PP", "Beta", "SE", "Z", "trait")
  trait_gr <- makeGRangesFromDataFrame(tab[,c(1,2,3,5,8,9)], keep.extra.columns = TRUE)
  ov1 <- findOverlaps(trait_gr, peak_gr)
  terminal_df <- data.frame(trait_gr[queryHits(ov1)], round(cpm[subjectHits(ov1),], 1))
  terminal_df[,c("seqnames", "start", "end", "PP", "trait", "HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8)))]
}

terminal_df <- lapply(c("HCT", "HGB", "MCH", "MCHC", "MEAN_RETIC_VOL", "RETIC_COUNT", "MCV", "RBC_COUNT"), import_to_df) %>%
  rbindlist() %>% data.frame() -> all_traits_snps

all_traits_snps %>%
  group_by(seqnames, start, end) %>%
  top_n(n = 1, wt = PP) %>% data.frame()  -> uniqueSNP

ATAC.cpm.log2.all <- uniqueSNP[, c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8)))]
ATAC.cpm.log2.all.mm.RBC <- ATAC.cpm.log2.all <- t(apply(ATAC.cpm.log2.all, 1, function(x)(x-min(x))/(max(x)-min(x))))

km <- kmeans(ATAC.cpm.log2.all.mm.RBC, centers = 9, nstart = 10000)
vals <- c("K4", "K6", "K7", "K9", "K5", "K1", "K8", "K3", "K2")
km.cluster <- factor(vals[as.numeric(km$cluster)], levels = sort(vals))

pdf(file="../../plots/ery_cluster_peaks_gwas.pdf", width = 3, height = 3)  
par(cex.main=0.8,mar=c(1,1,1,1))
hm <- Heatmap(ATAC.cpm.log2.all.mm.RBC, col=as.character(jdb_palette("solar_rojos",type="continuous")),
        cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 0),
        column_names_gp = gpar(fontsize = 6),
        split = km.cluster, show_heatmap_legend = FALSE,
        name = "Peak\nAccessibility")
hm
dev.off()

findBest <- cluster::clusGap(ATAC.cpm.log2.all.mm.RBC, FUN = kmeans, K.max = 15, B = 10, nstart = 1000)
qplot(1:15,findBest$Tab[,3])

uniqueSNP$Kcluster <- as.character(km.cluster)
# write.table(uniqueSNP, file = "../../processed/Kmeans_GWAS_lineage.tsv", sep = "\t", quote = FALSE, 
#             col.names = TRUE, row.names = FALSE)

# Make final scree plot
df <- data.frame(K = 1:15, Gap = findBest$Tab[,3])

p1 <- ggplot(df, aes(x = K, y = Gap)) + geom_point() +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "K-means", y = "Gap Statistic") +
  geom_vline(xintercept = 9, linetype = 2)

cowplot::ggsave(p1, file = "../../plots/gwas_peaks_kmeans_scree.pdf", width = 2.2, height =2)
