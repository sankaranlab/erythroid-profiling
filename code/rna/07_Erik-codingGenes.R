library(data.table)
library(dplyr)
library(BuenColors)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

set.seed(14651)
"%ni%" <- Negate("%in%")

lapply(list.files("../../processed/all_DEseq2/"), function(x){
  df <- read.table(paste0("../../processed/all_DEseq2/", x), header = TRUE)
  df$timepoint <- stringr::str_split_fixed(x, "_", 3)[1,2]
  df <- df[df$baseMean > 0,]
  df[,c("gene", "log2FoldChange", "timepoint")]
}) %>% rbindlist() %>% data.frame() -> allFC


lapply(list.files("../../downloads/RNA_DESeq2", full.names = TRUE), function(x){
  data.frame(g = read.table(x, header = TRUE, stringsAsFactors = FALSE)[,1])
}) %>% rbindlist() %>% data.frame() -> diffGeneDf

codingGenes <- read.table("../../data/erythroid_codinggenes_PP10.txt", stringsAsFactors = FALSE, header = FALSE)[,1]
go <- allFC[as.character(allFC$gene) %in% codingGenes, ]

square <- dcast(go, gene ~ timepoint, value.var = "log2FoldChange")
rownames(square) <- square[,1]
square <- data.matrix(square[,2:8])

km <- kmeans(square, centers = 5, nstart = 10000)
vals <- c("K4", "K5", "K1", "K3", "K2")
km.cluster <- factor(vals[as.numeric(km$cluster)], levels = sort(vals))

pdf(file="../../plots/ery_coding_peaks_gwas.pdf", width = 5, height = 3)  
par(cex.main=0.8,mar=c(1,1,1,1))
hm <- Heatmap(square, col=as.character(jdb_palette("brewer_spectra",type="continuous")),
        cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 0),
        column_names_gp = gpar(fontsize = 6),
        split = km.cluster, show_heatmap_legend = TRUE,
        name = "FC over P1")
hm
dev.off()


#codingGenes[codingGenes %ni% unique(go$gene)]
