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

# Intersect with GATA1 variants
x <- read.table("../../data/mendelianVariants.tsv", header = TRUE)
gata1_bp <- read.table("../../data/mendelianVariants_GATA1.tsv", header = FALSE)[,1]
x %>% filter(BP %in% gata1_bp) %>%
  mutate(start = BP, end = BP) %>% makeGRangesFromDataFrame() -> gata1mutgr

sh <- unique(subjectHits(findOverlaps(gata1mutgr, peak_gr)))

# Normalize
cpm <- sweep(raw, 2, colSums(raw), FUN="/") * 1000000
ATAC.cpm.log2.all <- log2(cpm[, c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8)))] +1)
ATAC.cpm.log2.all.mm.RBC <- t(apply(ATAC.cpm.log2.all[sh,], 1, function(x)(x-min(x))/(max(x)-min(x))))

pdf(file="../../plots/gata1_noncoding_peaks.pdf", width = 2, height = 2)  
par(cex.main=0.8,mar=c(1,1,1,1))
hm <- Heatmap(ATAC.cpm.log2.all.mm.RBC, col=as.character(jdb_palette("solar_rojos",type="continuous")),
        cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 0),
        column_names_gp = gpar(fontsize = 6),
        show_heatmap_legend = FALSE,
        name = "Peak\nAccessibility")
hm
dev.off()
