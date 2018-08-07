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
kmeans_df <- fread("../../processed/RNAseq-Kmeans7-clusterID-allpops.tsv")
outdir <- "../../data/FUMA/All_DiffExpressed/allpops/"
kmeans_df[kmeans_df$cluster == "K1",]$gene %>% write.table(.,paste0(outdir,"allDiff.k1.txt"),sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
kmeans_df[kmeans_df$cluster == "K2",]$gene %>% write.table(.,paste0(outdir,"allDiff.k2.txt"),sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
kmeans_df[kmeans_df$cluster == "K3",]$gene %>% write.table(.,paste0(outdir,"allDiff.k3.txt"),sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
kmeans_df[kmeans_df$cluster == "K4",]$gene %>% write.table(.,paste0(outdir,"allDiff.k4.txt"),sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
kmeans_df[kmeans_df$cluster == "K5",]$gene %>% write.table(.,paste0(outdir,"allDiff.k5.txt"),sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
kmeans_df[kmeans_df$cluster == "K6",]$gene %>% write.table(.,paste0(outdir,"allDiff.k6.txt"),sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
kmeans_df[kmeans_df$cluster == "K7",]$gene %>% write.table(.,paste0(outdir,"allDiff.k7.txt"),sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

symbol_keep %>% write.table(.,"../../data/FUMA/All_DiffExpressed/gene_universe.txt",sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Plot enrichments --------------------------------------------------------

plot_fuma <- function(toplot){
  toplot$GeneSet <- gsub("GO_","",toplot$GeneSet)
  toplot$GeneSet <- gsub("_"," ",toplot$GeneSet) %>% tolower
  toplot$GeneSet <- gsub("dna","DNA",toplot$GeneSet)
  toplot$GeneSet <- gsub("rna","RNA",toplot$GeneSet)
  toplot$GeneSet <- paste(toupper(substr(toplot$GeneSet, 1, 1)), substr(toplot$GeneSet, 2, nchar(toplot$GeneSet)), sep="")
  toplot$GeneSet <- factor(toplot$GeneSet,levels=toplot$GeneSet)
  
  ggplot(toplot,aes(y=-1*log(adjP),x=GeneSet,label=GeneSet)) + 
    geom_bar(stat="identity",fill=jdb_palette("solar_extra")[4],position = position_stack(reverse = TRUE)) + 
    coord_flip() +
    theme_bw() +
    pretty_plot(fontsize = 12) +
    L_border()+
    scale_y_continuous(expand=c(0,0))+
    scale_x_discrete(limits = rev(levels(toplot$GeneSet)))+
    theme(legend.position="none",
                axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    labs(x="",y="-log(FDR)") +
    geom_text(color="black",
              position = position_stack(vjust = 0.02),hjust=0,size=4)
}


indir <- "../../data/FUMA/All_DiffExpressed/allpops/"
outdir<-"../../tables/FUMA_allpops/"
top_n <- 100
for (i in 1:7){
  GS_k1 <- fread(paste0(indir,"GS_k",i,".txt")) %>%
    dplyr::select(-link) %>% filter(Category == "GO_bp") %>% arrange(adjP) %>% head(top_n) %>%
    dplyr::select(GeneSet,adjP,N_genes,N_overlap)-> toplot
  
  #Fix labels
  toplot$GeneSet <- gsub("GO_","",toplot$GeneSet)
  toplot$GeneSet <- gsub("_"," ",toplot$GeneSet) %>% tolower
  toplot <- toplot[toplot$GeneSet %in% goterms,]
  print(nrow(toplot)/nrow(GS_k1))
  toplot$GO_accession <- names(goterms[goterms %in% toplot$GeneSet])
  toplot$GeneSet <- gsub("dna","DNA",toplot$GeneSet)
  toplot$GeneSet <- gsub("rna","RNA",toplot$GeneSet)
  toplot$GeneSet <- paste(toupper(substr(toplot$GeneSet, 1, 1)), substr(toplot$GeneSet, 2, nchar(toplot$GeneSet)), sep="")
  toplot$GeneSet <- factor(toplot$GeneSet,levels=toplot$GeneSet)
  
  #write.table(toplot,file=paste0(outdir,"GO_enrichments_allpops_k",i,".tsv"),sep = "\t", quote = FALSE, col.names = T, row.names = FALSE)
  #toplot %>% dplyr::select(GO_accession,N_overlap,N_genes) %>% write.table(.,file=paste0(outdir,"totrim.k",i,".csv"),sep = ",", quote = FALSE, col.names = F, row.names = FALSE)
  if (i != 5){
    trimmed <- fread(paste0("../../tables/FUMA_allpops/k",i,"_trimmed.csv"))
    toplot <- toplot[trimmed$`Soft Trimmed?` != "YES",] 
  }
  plot_fuma(toplot[1:5,]) %>% ggsave(.,filename = paste0("../../plots/gene_set_enrichments/k",i,"_fuma_trimmed.pdf"), width =3.2, height = 2)
}

top_n=5
outdir <- "../../data/FUMA/All_DiffExpressed/allpops/"
GS_k1 <- fread(paste0(outdir,"GS_k1.txt")) %>%
  dplyr::select(-link) %>% filter(Category == "GO_bp") %>% arrange(adjP)
GS_k1 %>% head(top_n) %>% 
  plot_fuma() %>% ggsave(.,filename = "../../plots/gene_set_enrichments/k1_fuma.pdf", width =8, height = 4.5)

GS_k2 <- fread("../../data/FUMA/All_DiffExpressed/GS_k2.txt") %>%
  dplyr::select(-link) %>% filter(Category == "GO_bp") %>% arrange(adjP)
GS_k2 %>% head(top_n) %>% plot_fuma() %>%
  ggsave(.,filename = "../../plots/gene_set_enrichments/k2_fuma.pdf", width =8, height = 4.5)

GS_k3 <- fread("../../data/FUMA/All_DiffExpressed/GS_k3.txt") %>%
  dplyr::select(-link) %>% filter(Category == "GO_bp") %>% arrange(adjP)
GS_k3 %>% head(top_n) %>% plot_fuma() %>% 
  ggsave(.,filename = "../../plots/gene_set_enrichments/k3_fuma.pdf", width =8, height = 4.5)

GS_k4 <- fread("../../data/FUMA/All_DiffExpressed/GS_k4.txt") %>%
  dplyr::select(-link) %>% filter(Category == "GO_bp") %>% arrange(adjP)
GS_k4 %>% head(top_n) %>% plot_fuma() %>%
  ggsave(.,filename = "../../plots/gene_set_enrichments/k4_fuma.pdf", width =8, height = 4.5)

GS_k5 <- fread("../../data/FUMA/All_DiffExpressed/GS_k5.txt") %>%
  dplyr::select(-link) %>% filter(Category == "GO_bp") %>% arrange(adjP)
GS_k5 %>% head(top_n) %>% plot_fuma() %>% 
  ggsave(.,filename = "../../plots/gene_set_enrichments/k5_fuma.pdf", width =8, height = 4.5)
