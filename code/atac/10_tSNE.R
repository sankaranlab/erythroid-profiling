library(dplyr)
library(GenomicRanges)
library(BiocParallel)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(BuenColors)
library(Rtsne)
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
atac <- log2(cpm[, c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8)))] +1)
atac.norm <- t(apply(atac, 1, function(x)(x-min(x))/(max(x)-min(x))))

tsneO <- Rtsne(atac, check_duplicates = FALSE, perplexity = 10, pca = FALSE)
