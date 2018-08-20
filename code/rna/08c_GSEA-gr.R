library(liger)
library(data.table)
library(dplyr)
library(BuenColors)

mmtab <- read.table("../../data/mouseTargetGenes/nr3c1-1050362745114690_add6.tsv", stringsAsFactors = FALSE, header = TRUE)
klf1_mm <- mmtab %>% filter(adj.P.Val < 0.05 & logFC > 0) %>% pull(gene.symbols)

# mouse2human transition
lo_dt <- data.frame(fread("../../data/mouseTargetGenes/HOM_MouseHumanSequence.rpt", sep = "\t"))[,c(1,4)]

# Filter out non-pair concordances
lo_dt <- lo_dt[as.character(lo_dt$HomoloGene.ID) %in% attr(table(lo_dt$HomoloGene.ID), "names")[table(lo_dt$HomoloGene.ID) == 2], ]
stopifnot(all(lo_dt$HomoloGene.ID[c(TRUE,FALSE)] == lo_dt$HomoloGene.ID[c(FALSE,TRUE)]))
transdf <- data.frame(mouse = lo_dt[["Symbol"]][c(TRUE,FALSE)], human = lo_dt[["Symbol"]][c(FALSE,TRUE)], stringsAsFactors = FALSE)

mouse2human <- transdf[,"human"]; names(mouse2human) <- transdf[,"mouse"]

klf1_human <- unname(mouse2human[klf1_mm])[!is.na(mouse2human[klf1_mm])]

lapply(paste0("P", as.character(2:8)), function(pop){
  x <- data.frame(fread(paste0("../../processed/all_DEseq2/P1_",pop,"_deseq2full.tsv"))) %>% arrange(log2FoldChange)
  vals <- -1*x$log2FoldChange; names(vals) <- as.character(x$gene)
  data.frame(population = pop, edge =  bulk.gsea(values=vals, list(TARGETS= klf1_human), mc.cores=2)$edge)
}) %>% rbindlist() %>% data.frame() -> edgeDF

P1 <-  ggplot(edgeDF, aes(x = population, y = edge, group = 1)) + geom_line() + 
  geom_point() + geom_line() + pretty_plot(fontsize = 8) +
  L_border() + labs(x = "Differential expression population",
                    y = "pSTAT5 target genes\nGSEA edge statistic")

cowplot::ggsave(P1, file = "../../plots/GSEA-GR-targets.pdf", width = 2.4, height =2)

