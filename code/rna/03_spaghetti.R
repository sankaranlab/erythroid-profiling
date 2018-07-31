library(data.table)
library(dplyr)
library(BuenColors)

lapply(list.files("../../processed/all_DEseq2/"), function(x){
  df <- read.table(paste0("../../processed/all_DEseq2/", x), header = TRUE)
  df$timepoint <- stringr::str_split_fixed(x, "_", 3)[1,2]
  df <- df[df$baseMean > 500,]
  df[,c("gene", "log2FoldChange", "timepoint")]
}) %>% rbindlist() %>% data.frame() -> allFC


lapply(list.files("../../downloads/RNA_DESeq2", full.names = TRUE), function(x){
  data.frame(g = read.table(x, header = TRUE, stringsAsFactors = FALSE)[,1])
}) %>% rbindlist() %>% data.frame() -> diffGeneDf
differentialGenes <- unique(as.character(diffGeneDf[,1]))

go <- allFC[as.character(allFC$gene) %in% differentialGenes, ]

makePlot <- function(geneP){
  go <- rbind(go, data.frame(gene = unique(as.character(go$gene)), log2FoldChange = 0, timepoint = "P1"))
  go$interesting <- ifelse(go$gene == geneP, geneP, "zznone")
  go <- arrange(go, desc(interesting))
  go$gene <- factor(as.character(go$gene), levels = unique(as.character(go$gene)))
  
  colorVec <- c("firebrick", "grey"); names(colorVec) <- c(geneP,  "zznone")
  p1 <- ggplot(go, aes(x = timepoint, y = log2FoldChange*-1, group = gene, color = interesting)) +
    pretty_plot(fontsize = 8) +
    geom_line(size = 0.25) + scale_color_manual(values = colorVec) + L_border() +
    labs(x = "", y = "log2FC over P1") + theme(legend.position = "none")
  
  #cowplot::ggsave(p1, file = paste0("../../plots/spaghetti-",geneP,".pdf"), width = 3, height = 2)
  return(p1)
}
makePlot("TMCC2")
makePlot("GATA1")
makePlot("IRF1")
makePlot("EIF2AK1")
makePlot("HBB")


