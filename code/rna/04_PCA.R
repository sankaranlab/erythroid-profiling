library(dplyr)
library(irlba)
library(annotables)
library(BuenColors)
set.seed(14651)

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")

# Filter sex chromosomes
mdf <- merge(read.table("../../data/genes.tsv", header = FALSE), grch37, by.x = "V1", by.y = "ensgene")
symbol_keep <- unique(mdf[mdf$chr %in% c(as.character(1:22)),"V2"])

# Import raw expression values
raw <- read.table("../../processed/RNAseq_rawGeneCounts.tsv", header = TRUE)
raw <- raw[raw[,29] %in% symbol_keep,]
RNA.counts <- raw[,1:28]
meta <- stringr::str_split_fixed(colnames(RNA.counts), "_", 4)

# PCA setup
RNA.counts <- data.matrix(RNA.counts)
cpm <- sweep(RNA.counts, 2, colSums(RNA.counts), FUN="/") * 1000000
log2cpm <- log2(cpm + 1)
print(dim(log2cpm))

log2cpm_z <- sapply(1:dim(log2cpm)[1],function(i){
  x <- data.matrix(log2cpm)[i,]
  (x - mean(x))/sd(x)
}) %>% t()
log2cpm_z <- log2cpm_z[rowSums(is.finite(log2cpm_z)) == 28,]
pr_rna <- prcomp_irlba(log2cpm, n = 3, center = TRUE)
PCA_RNA <- data.frame(pr_rna$rotation, Population = meta[,3])

p1 <- ggplot(PCA_RNA, aes(x = PC1, y = PC2, color = Population)) + geom_point(size = 0.8) +
  scale_color_manual(values = eryth_color_maps) + pretty_plot(fontsize = 8) +
  theme( axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
   theme( axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + L_border() + labs(x = "RNA-seq PC1", y = "RNA-seq PC2") + 
  guides(color=guide_legend(ncol=2)) +theme(legend.position = "none")
  #theme(legend.position = c(0.7,0.25))

cowplot::ggsave(p1, file = "../../plots/RNA_PCA.pdf", width = 1.5, height = 1.5)
