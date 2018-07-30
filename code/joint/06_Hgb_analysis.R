library(data.table)
library(GenomicRanges)
library(dplyr)
library(diffloop)
library(qvalue)
        
# Look at HGB atac correlation peaks --------------------------------------

# Read output 
erycor.df <- fread(paste0("zcat < ","../../processed/peakGeneCorrelations_1Mb_all.tsv.gz"))
colnames(erycor.df) <- c("seqnames","start","end","gene","rho","p")
erycor.df$qvalue <- qvalue(erycor.df$p)$qvalues
erycor.df <- erycor.df %>%
  dplyr::filter(qvalue < 0.001)

# Calculate median number of unique peaks correlated with expression of each gene
temp <- erycor.df %>% group_by(gene) %>% summarize(n())
median(temp$`n()`)
plot(density(temp$`n()`))

# Calculate number of peaks associated with hemoglobin gene expression
# hgb_genes <- c("HBA1","HBA2","HBB","HBD","HBE1",
#                "HBG1","HBG2","HBM","HBQ1","HBZ")
hgb_genes <- c("HBA1","HBB")
temp2 <- erycor.df[erycor.df$gene %in% hgb_genes,] %>% group_by(gene) %>% summarize(n()) 
median(temp2$`n()`)

# Heatmap of Hgb-correlated peaks across heme states 

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

hgb_peaks <- erycor.df[erycor.df$gene %in% hgb_genes,] %>% makeGRangesFromDataFrame
#hgb_peaks <- pg.df[pg.df$gene %in% hgb_genes,] %>% makeGRangesFromDataFrame

ov1 <- findOverlaps(hgb_peaks, peak_gr)
hgb_countsonly <- round(cpm[subjectHits(ov1),], 1)
hgb_countsonly <- cpm[sample(nrow(cpm),size=1000,replace=FALSE),]

ATAC.cpm.log2.all.mm.Hgb <- t(apply(hgb_countsonly, 1, function(x)(x-min(x))/(max(x)-min(x))))

# Heatmap(ATAC.cpm.log2.all.mm.Hgb, col=as.character(jdb_palette("solar_rojos",type="continuous")),
#         cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
#         row_names_gp = gpar(fontsize = 0),
#         column_names_gp = gpar(fontsize = 6),
#         split = km.cluster, show_heatmap_legend = FALSE,
#         name = "Peak\nAccessibility")

Heatmap(ATAC.cpm.log2.all.mm.Hgb,cluster_columns=F,show_row_names = F)
