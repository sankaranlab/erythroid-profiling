library(data.table)
library(dplyr)
library(BuenColors)
library(circlize)
library(ComplexHeatmap)

# Import RNA-seq
rna <- data.frame(fread("../../processed/RNAseq_rawGeneCounts.tsv", header = TRUE))[,c(1:28)]
meta_RNA <- stringr::str_split_fixed(colnames(rna), "_", 4)
colnames(rna) <- paste0(meta_RNA[,3], "-", meta_RNA[,4]) 

# Import ATAC
atac <- data.frame(fread("../../data/ATAC_data/ery_only.counts.tsv", header = TRUE))[,c(1:28)]
meta_ATAC <- stringr::str_split_fixed(gsub("_ATAC.proatac", "", colnames(atac)), "_", 4)
colnames(atac) <- paste0(meta_ATAC[,3], "-", meta_ATAC[,4]) 

cor_rna <- cor(log2(rna + 1))

# Make plot
pdf(file="../../plots/ourData-RNAcor.pdf", width = 4, height = 3)  
par(cex.main=0.8,mar=c(1,1,1,1))
hm <- Heatmap(cor_rna[sort(rownames(cor_rna)), sort(rownames(cor_rna))], col=as.character(jdb_palette("brewer_spectra",type="continuous")),
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        show_heatmap_legend = TRUE,
        name = "Cor")
hm
dev.off()

cor_atac <- cor(atac)

# Make plot
pdf(file="../../plots/ourData-ATACcor.pdf", width = 4, height = 3)  
par(cex.main=0.8,mar=c(1,1,1,1))
hm <- Heatmap(cor_atac[sort(rownames(cor_atac)), sort(rownames(cor_atac))], col=as.character(jdb_palette("brewer_spectra",type="continuous")),
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        show_heatmap_legend = TRUE,
        name = "Cor")
hm
dev.off()


