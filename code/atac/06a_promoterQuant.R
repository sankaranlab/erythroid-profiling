library(dplyr)
library(GenomicRanges)
library(data.table)

# Import promoter annotations and gwas clusters
promoter <- diffloop::bedToGRanges("../../data/annotations/Promoter_UCSC.fixed.bed")
gwas <- read.table("../../processed/Kmeans_GWAS_lineage.tsv", header = TRUE)
gwas_gr <- makeGRangesFromDataFrame(gwas, keep.extra.columns = TRUE)

# Find overlaps and quant
ov <- findOverlaps(gwas_gr, promoter)
k_cluster_df <- data.frame(cluster = gwas$Kcluster, inPromoter = as.numeric(1:length(gwas_gr) %in% queryHits(ov)))
k_cluster_df %>% group_by(cluster) %>% summarize(propPromoter = mean(inPromoter))

