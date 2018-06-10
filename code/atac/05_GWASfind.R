library(dplyr)
library(GenomicRanges)
library(BiocParallel)
library(data.table)
register(MulticoreParam(2))

# Import
bed <- data.frame(fread("../../data/ATAC_data/ery_only.bed", header = FALSE))
colnames(bed) <- c("chr", "start", "end")
peak_gr <- makeGRangesFromDataFrame(bed)

# Import raw expression values
raw <- data.frame(fread("../../data/ATAC_data/ery_only.counts.tsv", header = TRUE))
ATAC.counts <- data.matrix(raw[,1:28])
meta <- stringr::str_split_fixed(colnames(ATAC.counts), "_", 4)

total <- sapply(unique(meta[,3]), function(pop){
  rowSums(ATAC.counts[,which(meta[,3] == pop)])
})

# Normalize
cpm <- sweep(total, 2, colSums(total), FUN="/") * 1000000

# Find overlaps with SNPS
import_to_GR <- function(file){
  tab <- read.table(file)
  tab <- tab[tab$V5 > 0.50,]
  colnames(tab) <- c("chr", "start", "end", "region", "PP", "Beta", "SE", "Z")
  gr <- makeGRangesFromDataFrame(tab[,c(1,2,3,5,8)], keep.extra.columns = TRUE)
}

trait_gr <- unique(import_to_GR("../../data/FMsnps/MPV_PP001_betas.bed"))
ov1 <- findOverlaps(RBCcount_gr, peak_gr)

terminal_df <- data.frame(trait_gr[queryHits(ov1)], round(cpm[subjectHits(ov1),], 1))
terminal_df <- data.frame(trait_gr[queryHits(ov2)], round(cpm[subjectHits(ov2),], 1))
terminal_df

