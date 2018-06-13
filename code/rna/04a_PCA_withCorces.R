library(dplyr)
library(irlba)
library(annotables)
library(BuenColors)
set.seed(14651)

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")
all_color_maps <- c(ejc_color_maps[c("HSC", "MPP", "CMP")], "MEP" = "#FF81AF" ,eryth_color_maps)

# Filter sex chromosomes
mdf <- merge(read.table("../../data/genes.tsv", header = FALSE), grch37, by.x = "V1", by.y = "ensgene")
symbol_keep <- unique(mdf[mdf$chr %in% c(as.character(1:22)),"V2"])

# Import raw expression values
load("../../processed/serialized/RNAseq_wCorces.rda")

# PCA setup
RNA.counts <- data.matrix(RNA.counts.all[,c(1:44)])
RNA.counts <- RNA.counts[RNA.counts.all[,45] %in% symbol_keep,]
cpm <- sweep(RNA.counts, 2, colSums(RNA.counts), FUN="/") * 1000000
log2cpm <- log2(cpm + 1)
print(dim(log2cpm))

log2cpm_z <- sapply(1:dim(log2cpm)[1],function(i){
  x <- data.matrix(log2cpm)[i,]
  (x - mean(x))/sd(x)
}) %>% t()

# Verify that genes aren't NAs
log2cpm_z <- log2cpm_z[rowSums(is.finite(log2cpm_z)) == 44,]
pr_rna <- prcomp_irlba(log2cpm_z, n = 3, center = TRUE)
PCA_RNA <- data.frame(pr_rna$rotation, Population = celltype)

p1 <- ggplot(PCA_RNA, aes(x = PC1, y = PC2, color = Population)) + geom_point(size = 0.8) +
  scale_color_manual(values = all_color_maps) + pretty_plot(fontsize = 8) +
 L_border() + labs(x = "RNA-seq PC1", y = "RNA-seq PC2") + 
  guides(color=guide_legend(ncol=3)) +theme(legend.position = "right")

cowplot::ggsave(p1, file = "../../plots/RNA_PCA_wCorces.pdf", width = 3.8, height = 2)
