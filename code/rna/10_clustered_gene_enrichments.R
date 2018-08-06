library(tidyverse)
library(annotables)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(plotly)
library(BuenColors)
library(qvalue)
library(lettercase)
"%ni%" <- Negate("%in%")

setwd("/Users/erikbao/Documents/GitHub/ery-final/code/atac")

# Load clustered differentially expressed genes  -----------------------------------------------
kmeans_df <- fread("../../processed/RNAseq-Kmeans5-clusterID.tsv")

kmeans_df[kmeans_df$cluster == "K1",]$gene %>% write.table(.,"../../data/FUMA/All_DiffExpressed/allDiff.k1.txt",sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
kmeans_df[kmeans_df$cluster == "K2",]$gene %>% write.table(.,"../../data/FUMA/All_DiffExpressed/allDiff.k2.txt",sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
kmeans_df[kmeans_df$cluster == "K3",]$gene %>% write.table(.,"../../data/FUMA/All_DiffExpressed/allDiff.k3.txt",sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
kmeans_df[kmeans_df$cluster == "K4",]$gene %>% write.table(.,"../../data/FUMA/All_DiffExpressed/allDiff.k4.txt",sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
kmeans_df[kmeans_df$cluster == "K5",]$gene %>% write.table(.,"../../data/FUMA/All_DiffExpressed/allDiff.k5.txt",sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

symbol_keep %>% write.table(.,"../../data/FUMA/All_DiffExpressed/gene_universe.txt",sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
kmeans_df[kmeans_df$cluster == "K1",]$gene -> k1genes

library(topGO)
# Filter sex chromosomes
mdf <- merge(read.table("../../data/genes.tsv", header = FALSE), grch37, by.x = "V1", by.y = "ensgene")
symbol_keep <- unique(mdf[mdf$chr %in% c(as.character(1:22), "X"),"V2"])
geneList <- factor(as.integer(symbol_keep %in% k1genes))

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
out <- getBM(filters= "hgnc_symbol",attributes= c("ensembl_gene_id","hgnc_symbol"),
      values=symbol_keep,mart= mart)

merge(df,G_list,by.x="gene",by.y="ensembl_peptide_id")

names(geneList) <- symbol_keep

geneID2GO <- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)

# Plot enrichments --------------------------------------------------------

plot_fuma <- function(toplot){
  toplot$GeneSet <- gsub("GO_","",toplot$GeneSet)
  toplot$GeneSet <- gsub("_"," ",toplot$GeneSet) %>% tolower
  toplot$GeneSet <- paste(toupper(substr(toplot$GeneSet, 1, 1)), substr(toplot$GeneSet, 2, nchar(toplot$GeneSet)), sep="")
  toplot$GeneSet <- factor(toplot$GeneSet,levels=toplot$GeneSet)
  
  ggplot(toplot,aes(y=-1*log(adjP),x=GeneSet)) + 
    geom_bar(stat="identity",fill="firebrick",position = position_stack(reverse = TRUE)) + 
    coord_flip() +
    theme_bw() +
    pretty_plot(fontsize = 8) +
    L_border()+
    scale_y_continuous(expand=c(0,0))+
    scale_x_discrete(limits = rev(levels(toplot$GeneSet)))+
    theme(legend.position="none")+
    labs(x="",y="-log(FDR)")
}

GS_k1 <- fread("../../data/FUMA/All_DiffExpressed/GS_k1.txt") %>%
  dplyr::select(-link) %>% filter(Category == "GO_bp") %>% arrange(adjP)
GS_k1 %>% head(10) %>% 
  plot_fuma() %>% ggsave(.,filename = "../../plots/gene_set_enrichments/k1_fuma.pdf", width =8, height = 4.5)

GS_k2 <- fread("../../data/FUMA/All_DiffExpressed/GS_k2.txt") %>%
  dplyr::select(-link) %>% filter(Category == "GO_bp") %>% arrange(adjP)
GS_k2 %>% head(10) %>% plot_fuma() %>%
  ggsave(.,filename = "../../plots/gene_set_enrichments/k2_fuma.pdf", width =8, height = 4.5)

GS_k3 <- fread("../../data/FUMA/All_DiffExpressed/GS_k3.txt") %>%
  dplyr::select(-link) %>% filter(Category == "GO_bp") %>% arrange(adjP)
GS_k3 %>% head(10) %>% plot_fuma() %>% 
  ggsave(.,filename = "../../plots/gene_set_enrichments/k3_fuma.pdf", width =8, height = 4.5)

GS_k4 <- fread("../../data/FUMA/All_DiffExpressed/GS_k4.txt") %>%
  dplyr::select(-link) %>% filter(Category == "GO_bp") %>% arrange(adjP)
GS_k4 %>% head(10) %>% plot_fuma() %>%
  ggsave(.,filename = "../../plots/gene_set_enrichments/k4_fuma.pdf", width =8, height = 4.5)

GS_k5 <- fread("../../data/FUMA/All_DiffExpressed/GS_k5.txt") %>%
  dplyr::select(-link) %>% filter(Category == "GO_bp") %>% arrange(adjP)
GS_k5 %>% head(10) %>% plot_fuma() %>% 
  ggsave(.,filename = "../../plots/gene_set_enrichments/k5_fuma.pdf", width =8, height = 4.5)
