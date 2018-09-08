library(dplyr)
library(BuenColors)
library(data.table)

# Import raw expression values
raw <- read.table("../../processed/RNAseq_rawGeneCounts.tsv", header = TRUE)
RNA.counts <- raw[,1:28]
meta <- stringr::str_split_fixed(colnames(RNA.counts), "_", 4)
genes <- as.character(raw[,29])
cpm <- sweep(RNA.counts, 2, colSums(RNA.counts), FUN="/") * 1000000
log2cpm <- log2(cpm + 1)

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")


# Population mean and SE
pops <- paste0("P", as.character(1:8))
lapply(pops, function(i){
  n <- sum(meta[,3] == i)
  means <- rowMeans(log2cpm[,meta[,3] == i])
  ses <- sqrt(matrixStats::rowVars(data.matrix(log2cpm[,meta[,3] == i])))/sqrt(n)
  data.frame(Population = i, 
             Gene = genes,
             Mean = means, 
             SEM = ses)
}) %>% rbindlist() %>% data.frame() -> odf

makeGeneBarPlot <- function(gene){
  odf %>% filter(Gene == gene) -> popAttribute
  
  popAttribute$Population <- factor(as.character(popAttribute$Population), 
                                    levels = rev(as.character(popAttribute$Population)))
  p1 <- ggplot(popAttribute, aes(x=Population, y=Mean, fill = Population)) + 
    geom_bar(position=position_dodge(), stat="identity", color = "black", width = 0.7) +
    coord_flip() + scale_fill_manual(values = eryth_color_maps) + 
    geom_errorbar(aes(ymin=Mean, ymax=Mean+SEM),
                  width=.2, position=position_dodge(.9)) +
    pretty_plot(fontsize = 8) + L_border() +
    labs(x = "", y = "") + theme(legend.position = "none") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    scale_y_continuous(breaks = pretty(popAttribute$Mean, n = 2), expand = c(0, 0))
  cowplot::ggsave(p1, file = paste0("../../plots/rna-bargraph/", gene, ".pdf"), 
                  width = 0.5, height = 2.02)  
}
makeGeneBarPlot("MYB")
makeGeneBarPlot("VEGFA")
makeGeneBarPlot("CCND3")
makeGeneBarPlot("UROS")
makeGeneBarPlot("RHAG")

makeGeneBarPlot("KLF1")
makeGeneBarPlot("KIT")
makeGeneBarPlot("HBA1")
makeGeneBarPlot("SMIM1")


