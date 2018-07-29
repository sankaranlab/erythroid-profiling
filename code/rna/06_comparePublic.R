library(data.table)
library(dplyr)
library(BuenColors)
library(circlize)
library(ComplexHeatmap)

ours <- data.frame(fread("../../processed/RNAseq_rawGeneCounts.tsv", header = TRUE))
theirs <- data.frame(fread(paste0("zcat < ", "../../data/GSE107218_CBPB-hg19-counts.txt.gz"), header = TRUE))
order_pub <- c("CB.CB_CD34_1","CB.CB_CD34_2","CB.CB_CD34_3","PB.PB_CD34_1","PB.PB_CD34_2","PB.PB_CD34_3",
               "CB.CB_BFUE_1","CB.CB_BFUE_2","CB.CB_BFUE_3","PB.PB_BFU_1","PB.PB_BFU_2","PB.PB_BFU_3",
               "CB.CB_CFUE_1","CB.CB_CFUE_2","CB.CB_CFUE_3","PB.PB_CFU_1","PB.PB_CFU_2","PB.PB_CFU_3",
               "CBPRO1","CBPRO2","CBPRO3","PB.Pro_PB_1","PB.Pro_PB_2","PB.Pro_PB_3","CBEBASO1","CBEBASO2",
               "CBEBASO3","PB.Early_baso_PB_1","PB.Early_baso_PB_2","PB.Early_baso_PB_3","CBLB1","CBLB2",
               "CBLB3","R.Late_baso_PB_1","R.Late_baso_PB_2","R.Late_baso_PB_3","CBPOLY1","CBPOLY3","CBPOLY2",
               "R.Poly_PB_1","R.Poly_PB_2","R.Poly_PB_3","CBORTHO1","CBORTHO2","CBORTHO3","R.Orth_PB_1","R.Orth_PB_2","R.Orth_PB_3")

theirs <- theirs[,c(1,7:54)]
mdf  <- merge(ours, theirs, by.x = "genes", by.y = "Geneid")

ours <- data.matrix(mdf[,2:29])
two_name <- stringr::str_split_fixed(colnames(ours), "_", 4)[,c(3,4)]
colnames(ours) <- sapply(1:dim(two_name)[1], function(i) paste0(two_name[i,1], "-", two_name[i,2]))

theirs <- data.matrix(mdf[,30:77])
colnames(theirs) <-   gsub("20120327.", "", gsub("20120309.", "", gsub("20111212.", "", gsub("hisat", "",gsub(".sorted.bam", "", colnames(theirs))))))
cormat <- cor(log2(ours+1), log2(theirs+1), method = "pearson")

cormat_plot <- cormat[ sort(rownames(cormat)), order_pub]

pdf(file="../../plots/compare_public_RNAdata.pdf", width = 5, height = 4)  
par(cex.main=0.8,mar=c(1,1,1,1))
hm <- Heatmap(cormat_plot, col=as.character(jdb_palette("brewer_spectra",type="continuous")),
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        show_heatmap_legend = FALSE,
        name = "Pearson Correlation")
hm
dev.off()

