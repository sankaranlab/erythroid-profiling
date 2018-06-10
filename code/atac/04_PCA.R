library(dplyr)
library(irlba)
library(annotables)
library(BuenColors)
set.seed(14651)

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")


# Import raw expression values
raw <- read.table("../../data/ATAC_data/ery_only.counts.tsv", header = TRUE)
ATAC.counts <- data.matrix(raw[,1:28])
meta <- stringr::str_split_fixed(colnames(ATAC.counts), "_", 4)

# PCA setup
cpm <- sweep(ATAC.counts, 2, colSums(ATAC.counts), FUN="/") * 1000000
log2cpm <- log2(cpm + 1)
print(dim(log2cpm))

log2cpm_z <- sapply(1:dim(log2cpm)[1],function(i){
  x <- data.matrix(log2cpm)[i,]
  (x - mean(x))/sd(x)
}) %>% t()
log2cpm_z <- log2cpm_z[rowSums(is.finite(log2cpm_z)) == 28,]
pr_atac <- prcomp_irlba(log2cpm_z, n = 3, center = TRUE)
PCA_ATAC <- data.frame(pr_atac$rotation, Population = meta[,3])

p1 <- ggplot(PCA_ATAC, aes(x = PC1, y = PC2, color = Population)) + geom_point(size = 0.8) +
  scale_color_manual(values = eryth_color_maps) + pretty_plot(fontsize = 8) +
  theme( axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
   theme( axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + L_border() + labs(x = "ATAC-seq PC1", y = "ATAC-seq PC2") + 
  guides(color=guide_legend(ncol=2)) +theme(legend.position = "none")
  #theme(legend.position = c(0.7,0.25))

cowplot::ggsave(p1, file = "../../plots/ATAC_PCA.pdf", width = 1.5, height = 1.5)
