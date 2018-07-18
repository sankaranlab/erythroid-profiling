library(data.table)
library(GenomicRanges)
library(dplyr)
library(diffloop)

# ---------------------------------------------
# Compute per-gene correlations with ATAC peaks
# ---------------------------------------------

# Import genomic coordinates
peaks_g <- makeGRangesFromDataFrame(fread("../../data/corces/panHeme.bed"),
                                    seqnames.field = "V1", start.field = "V2", end.field = "V3")
tss_g <- makeGRangesFromDataFrame(data.frame(fread("../../data/refGene_hg19_TSS.bed")), keep.extra.columns = TRUE,
                                  seqnames.field = "V1", start.field = "V2", end.field = "V2")

# Import atac
counts <-  data.matrix(data.frame(fread("../../data/corces/atac_combinedPopulations/panHeme.counts.tsv")))[,c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8)))]
ATAC.cpm <- round(sweep(counts, 2, colSums(counts), FUN="/") * 1000000, 1)
atac <- log2(ATAC.cpm+1)

# Import import RNA
load("../../processed/serialized/RNAseq_wCorces.rda")
RNA.counts <- data.matrix(RNA.counts.all[,c(1:44)])
RNAseq_counts <- sapply(c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8))), function(pop){
  rowSums(RNA.counts[,pop == celltype])
})
RNA.cpm <- round(sweep(RNAseq_counts, 2, colSums(RNAseq_counts), FUN="/") * 1000000, 1)
rna <- log2(RNA.cpm+1)
genes <- RNA.counts.all[,c(45)]

# And filter genes everywhere
keepGenes <- which(rowMeans(rna) > 1 & !is.na(match(genes, mcols(tss_g)$V5)))
length(keepGenes)
genes <- genes[keepGenes]
rna <- rna[keepGenes, ]
tss_g <- sort(padGRanges(tss_g[mcols(tss_g)[,3] %in% genes], pad = 1000000))

# Define a function to get pvalue from correlation given a sample size
rhoToP <- function(rho, n = 12){
  t <- abs(rho)*sqrt(n-2) / sqrt(1-rho^2)
  2*pt(t, n-1, lower=FALSE)
}

# Loop over genes to compute correlation, pvalue, etc. 
lappout <- lapply(genes, function(gene){
  atacidx <- unique(queryHits(findOverlaps(peaks_g, tss_g[which(mcols(tss_g)[,3] == gene)])))
  if(length(atacidx) > 1){
    rho <- cor((rna[which(genes == gene)[1],]), t((atac[atacidx,])))[1,]
    return(data.frame(data.frame(peaks_g[atacidx])[,c(1,2,3)], gene = gene, rho = rho, p = rhoToP(rho)))
  } else {
    return(NULL)
  }
})

# Make final output
df <- rbindlist(Filter(Negate(is.null), setNames(lappout,seq_along(lappout))))
df$rho <- round(df$rho, digits = 3)
df$p <- formatC(df$p, format = "e", digits = 2)
write.table(df, file = "../../processed/peakGeneCorrelations_1Mb_all.tsv", 
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)


