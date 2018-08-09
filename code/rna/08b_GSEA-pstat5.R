library(liger)
library(data.table)
library(dplyr)
library(BuenColors)

klf1_mm <- read.table("../../data/mouseTargetGenes/pSTAT5-tg.txt", stringsAsFactors = FALSE)[,1]

# mouse2human transition
lo_dt <- data.frame(fread("../../data/mouseTargetGenes/HOM_MouseHumanSequence.rpt", sep = "\t"))[,c(1,4)]

# Filter out non-pair concordances
lo_dt <- lo_dt[as.character(lo_dt$HomoloGene.ID) %in% attr(table(lo_dt$HomoloGene.ID), "names")[table(lo_dt$HomoloGene.ID) == 2], ]
stopifnot(all(lo_dt$HomoloGene.ID[c(TRUE,FALSE)] == lo_dt$HomoloGene.ID[c(FALSE,TRUE)]))
transdf <- data.frame(mouse = lo_dt[["Symbol"]][c(TRUE,FALSE)], human = lo_dt[["Symbol"]][c(FALSE,TRUE)], stringsAsFactors = FALSE)

mouse2human <- transdf[,"human"]; names(mouse2human) <- transdf[,"mouse"]

klf1_human <- unname(mouse2human[klf1_mm])[!is.na(mouse2human[klf1_mm])]

x <- data.frame(fread("../../processed/all_DEseq2/P1_P8_deseq2full.tsv")) %>% arrange(log2FoldChange)
vals <- -1*x$log2FoldChange; names(vals) <- as.character(x$gene)
bulk.gsea(values=vals, list(KLF1_TARGETS= klf1_human), mc.cores=2) 


edgeDF <- data.frame(
  edge = c(-4.1,-5.9,-1.7,1.4,1.5,3.1,3.4),
  # sscore = c(1.779383, 1.733049,  1.911498, 2.219836, 2.101272, 2.021585, 1.726593),
  population = c("P2", "P3", "P4", "P5", "P6", "P7", "P8")
)

P1 <-  ggplot(edgeDF, aes(x = population, y = edge, group = 1)) + geom_line() + 
  geom_point() + geom_line() + pretty_plot(fontsize = 8) +
  L_border() + labs(x = "Differential expression population",
                    y = "pSTAT5 target genes\nGSEA edge statistic")

cowplot::ggsave(P1, file = "../../plots/GSEA-pSTAT5-targets.pdf", width = 2.4, height =2)

