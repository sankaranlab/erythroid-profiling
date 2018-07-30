library(reshape2)
library(dplyr)
library(BuenColors)
library(GenomicRanges)
library(data.table)

# Import loops 
cd34_loops <- readRDS("../../data/loops/cd34_loops.rds")
ery_loops <- readRDS("../../data/loops/ery_loops.rds")
ery_loops$chr <- paste0("chr",ery_loops$oeChr)
ery_loops$which <- "Ery"
cd34_loops$which <- "CD34"

# Combined DF
df <- data.frame(
  chr = c(ery_loops$chr, cd34_loops$chr),
  start = c(ery_loops$oeStart, cd34_loops$start),
  end = c(ery_loops$oeEnd, cd34_loops$end),
  gene = c(ery_loops$baitName, cd34_loops$gene),
  which = c(ery_loops$which, cd34_loops$which)
)

gr <- makeGRangesFromDataFrame(df,keep.extra.columns = TRUE)

# Import ATAC & normalize
peaksdf <- fread("../../data/corces/panHeme.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(data.frame(fread("../../data/corces/atac_combinedPopulations/panHeme.counts.tsv")))[,c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8)))]
ATAC.cpm <- round(sweep(counts, 2, colSums(counts), FUN="/") * 1000000, 1)
ATAC.cpm.log2 <- log2(ATAC.cpm+1)

# Handle RNA-seq
load("../../processed/serialized/RNAseq_wCorces.rda")
RNA.counts <- data.matrix(RNA.counts.all[,c(1:44)])
RNAseq_counts <- sapply(c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8))), function(pop){
  rowSums(RNA.counts[,pop == celltype])
})
RNA.cpm <- round(sweep(RNAseq_counts, 2, colSums(RNAseq_counts), FUN="/") * 1000000, 1)
RNA.cpm.log2 <- log2(RNA.cpm+1)

# Exclude P1
RNA.cpm.log2 <- RNA.cpm.log2[,!(colnames(RNA.cpm.log2) == "P1")]
ATAC.cpm.log2 <- ATAC.cpm.log2[,!(colnames(ATAC.cpm.log2) == "P1")]

# Find maximum populations
gene_max <- data.frame(gene = RNA.counts.all$genes, maxPopRNA = colnames(RNA.cpm.log2)[max.col( RNA.cpm.log2)])
peak_max <- data.frame(peaks, maxPopATAC = colnames(ATAC.cpm.log2)[max.col(ATAC.cpm.log2)])

# Intersect the two via pchic loops
ov <- findOverlaps(peaks, gr)
gene_enh <- data.frame(peak_max[queryHits(ov),], mcols(gr)[subjectHits(ov),])

alltogether <- merge(gene_enh, gene_max)

alltogether %>% group_by(maxPopRNA, maxPopATAC) %>% summarise(count = n()) %>% data.frame() -> quants

ggplot(quants, aes(x = maxPopRNA, y = maxPopATAC, fill = count)) + geom_tile() + pretty_plot() + 
  scale_fill_gradientn(colors = jdb_palette("solar_rojos")) 
#  geom_abline(slope=-1, intercept = 9, linetype = 2)  + labs(x = "", y = "") + theme(legend.position = "none")


