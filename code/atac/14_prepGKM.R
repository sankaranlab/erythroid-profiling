library(BuenColors)
library(dplyr)
library(DESeq2)
library(BiocParallel)
library(data.table)
library(annotables)

register(MulticoreParam(2))

# Import raw expression values
ATAC.counts <- read.table("../../data/ATAC_data/combined_12.counts.tsv", header = TRUE)[,c(5:12)]
peaks <- data.frame(fread("../../data/ATAC_data/ery_only.bed"))

scaleRows <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
}

# get quantiles for Z scores
ATAC.counts.filt.Z <- scaleRows(ATAC.counts)
Qs <- apply(ATAC.counts.filt.Z,2,function(x) {quantile(x,0.90)})
Ps <- paste0("P", as.character(1:8))

lapply(1:8, function(i){
  write.table(peaks[which(ATAC.counts.filt.Z[,i] > Qs[i]),], 
              file = paste0("../../processed/gkmerpeaks/topPeaks-", Ps[i], ".bed"),
              row.names = FALSE, col.names = FALSE, quote = FALSE) 
})

