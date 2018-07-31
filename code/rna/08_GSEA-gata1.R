library(liger)
library(data.table)
library(dplyr)
library(BuenColors)

load("../../data/annotations/MSigDB_Hs_GS.RData")
ll <- ls(MSigDB.Hs.GS2Symbol)
geneSets <- sapply(ll, function(x) get(x, env=MSigDB.Hs.GS2Symbol))
gs <- geneSets["WELCH_GATA1_TARGETS"]

x <- data.frame(fread("../../processed/all_DEseq2/P1_P5_deseq2full.tsv")) %>% arrange(log2FoldChange)
vals <- -1*x$log2FoldChange; names(vals) <- as.character(x$gene)
gsea(values=vals, gs, mc.cores=2, plot = TRUE) 

x <- data.frame(fread("../../processed/all_DEseq2/P1_P7_deseq2full.tsv")) %>% arrange(log2FoldChange)
vals <- -1*x$log2FoldChange; names(vals) <- as.character(x$gene)
bulk.gsea(values=vals, gs, mc.cores=2) 


 names(geneSets)[grep("GATA", names(geneSets))]
 
 edgeDF <- data.frame(
   edge = c(0.8, 1.3, 1.2, 3.2, 3.2, 2.9, 2.2),
   sscore = c(1.779383, 1.733049,  1.911498, 2.219836, 2.101272, 2.021585, 1.726593),
   population = c("P2", "P3", "P4", "P5", "P6", "P7", "P8")
 )
 
 P1 <-  ggplot(edgeDF, aes(x = population, y = edge, group = 1)) + geom_line() + 
   geom_point() + geom_line() + pretty_plot(fontsize = 8) +
   L_border() + labs(x = "Differential expression population",
                     y = "GATA1 target genes\nGSEA edge statistic")
 
 cowplot::ggsave(P1, file = "../../plots/GSEA-GATA1-targets.pdf", width = 2.4, height =2)
 
 