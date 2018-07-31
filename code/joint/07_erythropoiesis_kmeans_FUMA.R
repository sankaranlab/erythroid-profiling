library(tidyverse)
library(annotables)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(plotly)
library(BuenColors)
library(qvalue)
"%ni%" <- Negate("%in%")

#setwd("/Users/erikbao/Documents/GitHub/ery-final/")

# Load kmeans GWAS variants -----------------------------------------------

# Fine-map variants PP>0.10
kmeans_df <- fread("../../processed/Kmeans_GWAS_lineage.tsv")
CS.gr <- readRDS("/Users/erikbao/Documents/GitHub/singlecell_bloodtraits/data/Finemap/UKBB_BC_v3_VEPannotations.rds")
CS.df <- CS.gr %>% as.data.frame() 
CS.df$merge <- paste(CS.df$seqnames,CS.df$start,sep=":")
kmeans_df$merge <- paste(kmeans_df$seqnames,kmeans_df$start,sep=":")
merged <- merge(kmeans_df,CS.df[,c("merge","var","Consequence","SYMBOL")],by="merge") %>% dplyr::select(-merge) %>% unique()

m1 <- merged[merged$Kcluster %in% c("K6","K7"),] %>% group_by(trait) %>% summarise(n())
m2 <- merged %>% group_by(trait) %>% summarise(n())


# PChIC Gene Targets ------------------------------------------------------
# Read in PCHIC 
pchic <- readRDS("../../..//singlecell_bloodtraits/data/pchic/pchic.rds")
pchic.gr <- GRanges(pchic)
grch38.pc <- grch38 %>%
  dplyr::filter(biotype == "protein_coding")
pchic.gr <- pchic.gr[pchic.gr$Gene %in% grch38.pc$symbol,]
# # Read in CD34 PCHIC
# pchic_cd34 <- readRDS("/Volumes/broad_sankaranlab/ebao/singlecell_bloodtraits/data/pchic/cd34_loops.rds")
# pchic_cd34.gr <- GRanges(pchic_cd34)
# pchic_cd34.gr <- pchic_cd34.gr[pchic_cd34.gr$gene %in% grch38.pc$symbol,]

erytraits <- c("HCT", "HGB", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL", "RBC_COUNT", "RETIC_COUNT")
pchic.gr.RBC <- pchic.gr[pchic.gr$CellType == "Ery" & pchic.gr$PP > 0.1 & pchic.gr$Trait %in% erytraits,]
pchic.df.RBC <- as.data.frame(pchic.gr.RBC)
pchic.df.RBC$tomerge <- paste(pchic.df.RBC$seqnames,pchic.df.RBC$variantPos,sep=":")
merged$tomerge <- paste(merged$seqnames,merged$start,sep=":")
pchic.merged <- left_join(merged,pchic.df.RBC[,c("tomerge","Gene","Value","rsid")],by="tomerge") %>% dplyr::select(-tomerge)

pchic.merged[pchic.merged$Kcluster %in% c("K6","K7") &
               pchic.merged$trait == "HGB",] %>% dplyr::select(PP,var,Gene) %>% distinct()


# ATAC-RNA correlations gene targets from larger hemeATAC --------------------------------------
# Read in enhancer gene correlations
pg.df <- fread(paste0("zcat < ", "../../../singlecell_bloodtraits/data/bulk/peakGeneCorrelation.tsv.gz"))
names(pg.df) <- c("chrom","j.start","j.end","gene","cor","pvalue")
pg.df$qvalue <- qvalue(pg.df$pvalue)$qvalues
pg.df <- pg.df %>%
  dplyr::filter(qvalue < 0.001)

setkey(setDT(merged), seqnames, start, end)
setkey(setDT(pg.df), chrom, j.start, j.end)
ATAC.cor <- foverlaps(merged,pg.df, nomatch = 0)  %>%
  deployer::filter(PP > 0.1) %>%
  dplyr::select(var, trait, PP, seqnames, start, end, gene, cor, qvalue,Kcluster)

ATAC.cor.RBC <- ATAC.cor %>%
  dplyr::filter(trait %in% c("HCT", "HGB", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL", "RBC_COUNT", "RETIC_COUNT")) %>%
  GRanges()

make_atac_genelist <- function(df,clust){
  filtered_cs <- df %>% as.data.frame() %>% subset(Kcluster %in% clust)
  return(filtered_cs$gene %>% na.omit %>%unique)
}

genelist1 <- make_atac_genelist(ATAC.cor.RBC,clust=paste0("K",1:4))
write.table(genelist, file = "data/FUMA/k1_4_pchic_genes.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

genelist <- make_atac_genelist(ATAC.cor.RBC,clust=paste0("K",8:9))
write.table(genelist, file = "data/FUMA/k8_9_atac_cor_genes.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

ATAC.cor.RBC %>% subset(gene=="GFI1B")


# ATAC-RNA correlations from erythroid ATAC -------------------------------

# Read dataset 
erycor.df <- fread(paste0("zcat < ","../../processed/peakGeneCorrelations_1Mb_all.tsv.gz"))
colnames(erycor.df) <- c("seqnames","start","end","gene","rho","p")
erycor.df$qvalue <- qvalue(erycor.df$p)$qvalues
erycor.df <- erycor.df %>%
  dplyr::filter(qvalue < 0.001)

names(erycor.df) <- c("chrom","j.start","j.end","gene","cor","pvalue","qvalue")
setkey(setDT(merged), seqnames, start, end)
setkey(setDT(erycor.df), chrom, j.start, j.end)
ATAC.cor <- data.table::foverlaps(merged,erycor.df, nomatch = 0)  %>%
  deployer::filter(PP > 0.1) %>%
  dplyr::select(var, trait, PP, seqnames, start, end, gene, cor, qvalue,Kcluster)

ATAC.cor.RBC <- ATAC.cor %>%
  dplyr::filter(trait %in% c("HCT", "HGB", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL", "RBC_COUNT", "RETIC_COUNT")) %>%
  GRanges()
# FUMA PChic gene to function --------------------------------------------------------
make_pchic_genelist <- function(df,clust){
  filtered_cs <- df %>% subset(Kcluster %in% clust)
  return(filtered_cs$Gene %>% na.omit %>%unique)
}

genelist2 <- make_pchic_genelist(pchic.merged,clust=paste0("K",1:4))
write.table(genelist, file = "data/FUMA/k1_4_pchic_genes.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

genelist <- make_pchic_genelist(pchic.merged,clust=paste0("K",5:7))
write.table(genelist, file = "data/FUMA/k5_7_pchic_genes.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

genelist <- make_pchic_genelist(pchic.merged,clust=paste0("K",8:9))
write.table(genelist, file = "data/FUMA/k8_9_pchic_genes.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

pchic.merged %>% subset(Kcluster %in% paste0("K",8:9)) %>% na.omit %>% filter(Gene =="JAK2")
pchic.merged %>% subset(Kcluster %in% paste0("K",8:9)) %>% na.omit %>% filter(Gene =="CALR")
pchic.merged %>% subset(Kcluster %in% paste0("K",5:7)) %>% na.omit %>% filter(Gene =="JAK2")


# FUMA coding gene to function --------------------------------------------------------------------

make_coding_genelist <- function(df,clust){
  coding_consequences <- c("missense_variant","synonymous_variant","frameshift_variant",
                           "splice_acceptor_variant","splice_donor_variant","splice_region_variant",
                           "inframe_insertion","stop_gained","stop_retained_variant",
                           "start_lost","stop_lost","coding_sequence_variant","incomplete_terminal_codon_variant")
  filtered_cs <- df %>% subset(Kcluster %in% clust) %>% filter(SYMBOL != "-") %>% subset(Consequence %in% coding_consequences) 
  return(filtered_cs$SYMBOL %>% na.omit %>%unique)
}

merged %>% subset(Kcluster %in% paste0("K",8:9))  %>% filter(SYMBOL != "-")
make_coding_genelist(merged,clust=paste0("K",1:4))
make_coding_genelist(merged,clust=paste0("K",5:7))
make_coding_genelist(merged,clust=paste0("K",8:9))


# FUMA with all target genes ----------------------------------------------

K1K4 <- dplyr::union_all(make_atac_genelist(ATAC.cor.RBC,clust=paste0("K",1:4)), 
              make_pchic_genelist(pchic.merged,clust=paste0("K",1:4)),
              make_coding_genelist(merged,clust=paste0("K",1:4)))
K5K7 <- dplyr::union_all(make_atac_genelist(ATAC.cor.RBC,clust=paste0("K",5:7)), 
                         make_pchic_genelist(pchic.merged,clust=paste0("K",5:7)),
                         make_coding_genelist(merged,clust=paste0("K",5:7)))
K8K9 <- dplyr::union_all(make_atac_genelist(ATAC.cor.RBC,clust=paste0("K",8:9)), 
                         make_pchic_genelist(pchic.merged,clust=paste0("K",8:9)),
                         make_coding_genelist(merged,clust=paste0("K",8:9)))
write.table(K8K9, file = "data/FUMA/k8_9_all_target_genes.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)


# MotifbreakR -------------------------------------------------------------
library(motifbreakR)
motifbreakr <- readRDS("../singlecell_bloodtraits/data/motifbreakR/alltraits.mbreaker.withPPs.rds")
motifbreakr[motifbreakr$SNP %in% "chr20:4162978:C:T",]

TMCC2_coordinates <- c(205197008,205242501)
flank <- 10000
CS.gr %>% subset(seqnames == "chr1" & 
                   start > (TMCC2_coordinates[1]-flank) & 
                   end < TMCC2_coordinates[2]+flank)
