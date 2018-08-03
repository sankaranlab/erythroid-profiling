library(tidyverse)
library(annotables)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(plotly)
library(BuenColors)
library(qvalue)
"%ni%" <- Negate("%in%")

setwd("/Users/erikbao/Documents/GitHub/ery-final/code/atac")

# Load kmeans GWAS variants -----------------------------------------------

# Fine-map variants PP>0.10
kmeans_df <- fread("../../processed/Kmeans_GWAS_lineage.tsv")
kmeans_df_nonunique <- fread("../../processed/Kmeans_GWAS_lineage_nonunique.tsv")
joined <- left_join(kmeans_df_nonunique,kmeans_df[,c("seqnames","start","Kcluster")],by=c("seqnames","start")) %>% distinct()

CS.gr <- readRDS("/Users/erikbao/Documents/GitHub/singlecell_bloodtraits/data/Finemap/UKBB_BC_v3_VEPannotations.rds")
CS.df <- CS.gr %>% as.data.frame() 
CS.df$merge <- paste(CS.df$seqnames,CS.df$start,sep=":")
joined$merge <- paste(joined$seqnames,joined$start,sep=":")

merged <- merge(joined,CS.df[,c("merge","var","Consequence","SYMBOL")],by="merge") %>% dplyr::select(-merge) %>% unique()

# Number of variants in each cluster/group of clusters
m1 <- merged[merged$Kcluster %in% c("K6","K7"),] %>% group_by(trait) %>% summarise(n())
m2 <- merged %>% group_by(trait) %>% summarise(n())

# Variants that are in peaks stronger in P2-P4
merged[merged$Kcluster %in% c("K6") & merged$trait %in% "HGB",]  %>% distinct(var)

# MotifbreakR -------------------------------------------------------------
library(motifbreakR)
motifbreakr <- readRDS("../../../singlecell_bloodtraits/data/motifbreakR/alltraits.mbreaker.withPPs.rds")
motifbreakr[motifbreakr$SNP %in% "chr20:4162978:C:T",]

TMCC2_coordinates <- c(205197008,205242501)
flank <- 10000
CS.gr %>% subset(seqnames == "chr1" & 
                   start > (TMCC2_coordinates[1]-flank) & 
                   end < TMCC2_coordinates[2]+flank)
