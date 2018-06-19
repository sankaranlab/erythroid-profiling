library(data.table)
library(dplyr)
library(BuenColors)
library(circlize)
library(ComplexHeatmap)

ours <- data.frame(fread("../../processed/RNAseq_rawGeneCounts.tsv", header = TRUE))
theirs <- data.frame(fread(paste0("zcat < ", "../../data/GSE107218_CBPB-hg19-counts.txt.gz"), header = TRUE))
theirs <- theirs[,c(1,7:54)]
mdf  <- merge(ours, theirs, by.x = "genes", by.y = "Geneid")

ours <- data.matrix(mdf[,2:29])
two_name <- stringr::str_split_fixed(colnames(ours), "_", 4)[,c(3,4)]
colnames(ours) <- sapply(1:dim(two_name)[1], function(i) paste0(two_name[i,1], "-", two_name[i,2]))

theirs <- data.matrix(mdf[,30:77])
colnames(theirs) <-   gsub("20120327.", "", gsub("20120309.", "", gsub("20111212.", "", gsub("hisat", "",gsub(".sorted.bam", "", colnames(theirs))))))
cormat <- cor(log2(ours+1), log2(theirs+1), method = "pearson")

pdf(file="../../plots/compare_public_RNAdata.pdf", width = 5, height = 4)  
par(cex.main=0.8,mar=c(1,1,1,1))
hm <- Heatmap(cormat, col=as.character(jdb_palette("brewer_spectra",type="continuous")),
        cluster_rows = TRUE, cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        show_heatmap_legend = FALSE,
        name = "Pearson Correlation")
hm
dev.off()

