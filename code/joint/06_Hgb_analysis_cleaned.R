# Libraries
library(BuenColors)
library(data.table)
library(GenomicRanges)
library(reshape2)
library(ComplexHeatmap)
library(matrixStats)
library(SummarizedExperiment)
library(Matrix)
library(preprocessCore)
library(tidyverse)
library(annotables)
library(qvalue)
library(chromVAR)
library(gchromVAR)
library(chromVARxx)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggextra)
library(ggrepel)
"%ni%" <- Negate("%in%")

# Identifying g-chromVAR trends -------------------------------------------

# Import and run run run 
peaksdf <- fread("../../data/ATAC_data/ery_only.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
ukbb <- importBedScore(peaks, list.files("../../data/FMsnps/", full.names = TRUE, pattern = "*.bed$"))
ukbb <- ukbb[keep,]

ATAC.counts <- data.matrix(data.frame(data.table::fread("../../data/ATAC_data/ery_only.counts.tsv")))
meta <- stringr::str_split_fixed(colnames(ATAC.counts), "_", 4)
pops.ordered <- c( paste0("P", as.character(1:8)))

sapply(pops.ordered, function(pop){
  rowSums(ATAC.counts[,meta[,3] == pop])
}) -> counts

if (TRUE){
  # Keep only good peaks in nth percentile of at least one cell type
  n = 0.8
  counts_normalized <- normalize.quantiles(as.matrix(counts)) # Quantile normalize distributions across 8 pops
  keep <- apply(counts_normalized,1,max) > mean(apply(counts_normalized,2,function(x) {quantile(x,n)})) # Compare max count in a peak with the mean nth percentile of all peaks per cell
  counts <- counts[keep,]
  
  # Subset to only good peaks
  peaks <- peaks[keep,]
}

# Create objects for g-chromVAR
cpm <- round(sweep(counts, 2, colSums(counts), FUN="/") * 1000000, 1)
SE <- SummarizedExperiment(assays = list(counts = counts),
                           rowData = peaks, 
                           colData = DataFrame(names = colnames(counts)))
SE <- filterPeaks(SE)
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)

# Take only peaks with a cumulative PP of x across all erythroid traits
erytraits <- c("HCT", "HGB", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL", "RBC_COUNT", "RETIC_COUNT")
keepPeaks <- rowSums(assays(ukbb)[["weights"]][,paste0(erytraits, "_PP001_betas")]) > 0.1
mega_df <- data.frame(
  rowRanges(SE)[keepPeaks],
  data.matrix(cpm[keepPeaks,]),
  data.matrix(assays(ukbb)[["weights"]][keepPeaks,paste0(erytraits, "_PP001_betas")]))
colnames(mega_df)[15:ncol(mega_df)] <- erytraits

mega_df$meanP2_4 <- apply(mega_df[,paste0("P",seq(2,4))],1,mean)
mega_df$meanP5_6 <- apply(mega_df[,paste0("P",seq(5,6))],1,mean)
mega_df$FC <- log2(mega_df$meanP5_6  / mega_df$meanP2_4)
#mega_df$FC <- log2(mega_df$P6  / mega_df$P4)


# Overlap with finemapped variants
CS.gr <- readRDS("../../data/UKBB_BC_v3_VEPannotations.rds")
CS.df <- CS.gr %>% as.data.frame() %>% as.data.table
setkey(setDT(CS.df), seqnames, start, end)

mega_df<-as.data.table(mega_df)
names(mega_df) <- c("chrom","j.start","j.end",names(mega_df)[4:ncol(mega_df)])
setkey(setDT(mega_df), chrom, j.start, j.end)


# HGB/HCT plots -----------------------------------------------------------
traits_to_plot <- c("HGB","HCT")
hgb_variants <- foverlaps(CS.df,mega_df, nomatch = 0)  %>%
  dplyr::filter(trait %in% traits_to_plot) %>% arrange(desc(PP)) %>% 
  dplyr::select(var,trait, seqnames, start,end,PP,meanP2_4,meanP5_6,FC) %>%
  distinct(var,.keep_all=TRUE)

# HGB/HCT Volcano plot
selected_vars <- c("6:43742626_T_C","6:43737805_A_C",
                   "4:55408875_A_T","4:55408999_T_C",
                   "1:3691528_A_G") #VEGFA, KIT, SMIM1
hgb_variants$highlight<- "F"
#hgb_variants[hgb_variants$var %in% selected_vars,"highlight"] <- "T"

p1 <-ggplot(hgb_variants, aes(x=FC, y=PP)) +
  geom_point(aes(color=highlight)) +
  scale_color_manual(values = jdb_palette("solar_extra")[c(1,7)]) +
  pretty_plot(fontsize = 10) + 
  L_border() + 
  geom_vline(xintercept = 0, linetype = 2) + 
  labs(x="log2FC",y="Posterior Probability") + 
  theme(legend.position="none")

# ggExtra::ggMarginal(p1, type = "density",margins="x")
#cowplot::ggsave(p1, file="../../plots/HGB_HCT_variants.pdf", width =3.5, height = 2.5)

# Weighted density plot
d1 <-ggplot(hgb_variants, aes(FC)) +
  geom_density(aes(weight=abs(FC)*PP/sum(abs(FC)*PP))) +
  pretty_plot(fontsize = 10)+
  L_border() +
  labs(x="")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.33)) +
  geom_vline(xintercept = 0, linetype = 2) 
#cowplot::ggsave(d1, file="../../plots/HGB_HCT_variants_weighteddensity.pdf", width =3.5, height = 1)


# MCV/MCH/MCHC/RETIC/RBC plots --------------------------------------------
traits_to_plot <- c("MCV","MCH","RBC_COUNT","MCHC")
mcv_variants <- foverlaps(CS.df,mega_df, nomatch = 0)  %>%
  dplyr::filter(trait %in% traits_to_plot) %>% arrange(desc(PP)) %>% 
  dplyr::select(var,trait, seqnames, start,end,PP,meanP2_4,meanP5_6,FC) %>%
  distinct(var,.keep_all=TRUE)

# Volcano plot
selected_vars <- c("19:12994618_C_T",
                   "6:41924998_C_T","6:41925159_G_A",
                   "16:170076_G_C") # KLF1, CCND3, HBA1
mcv_variants$highlight<- "F"
#mcv_variants[mcv_variants$var %in% selected_vars,"highlight"] <- "T"
p2 <-ggplot(mcv_variants, aes(x=FC, y=PP)) +
  geom_point(aes(color=highlight)) +
  scale_color_manual(values = jdb_palette("solar_extra")[c(1,7)]) +
  pretty_plot(fontsize = 10) + 
  L_border() + 
  theme(legend.position="none")+
  geom_vline(xintercept = 0, linetype = 2) + 
  labs(x="log2FC",y="Posterior Probability") 

#cowplot::ggsave(p2, file="../../plots/MCV_MCH_MCHC_RBC_variants.pdf", width =3.5, height = 2.5)

# Weighted density
d2 <- ggplot(mcv_variants, aes(FC)) +
  geom_density(aes(weight=abs(FC)*PP/sum(abs(FC)*PP))) +
  pretty_plot(fontsize = 10)+
  L_border() +
  labs(x="log2FC")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.33)) +
  geom_vline(xintercept = 0, linetype = 2) 
#cowplot::ggsave(d2, file="../../plots/MCV_MCH_MCHC_RBC_variants_weighteddensity.pdf", width =3.5, height = 1)


# Merge with ATAC-RNA correlations ----------------------------------------
pg.df <- fread(paste0("zcat < ", "/Users/erikbao/Documents/GitHub/singlecell_bloodtraits/data/bulk/peakGeneCorrelation.tsv.gz"))
names(pg.df) <- c("chrom","j.start","j.end","gene","cor","pvalue")
pg.df$qvalue <- qvalue(pg.df$pvalue)$qvalues
pg.df <- pg.df %>%
  dplyr::filter(qvalue < 0.001)
setkey(setDT(pg.df), chrom, j.start, j.end)

hgb_variants <- as.data.table(hgb_variants)
setkey(setDT(hgb_variants), seqnames, start, end)

hgb_atac_cor<- foverlaps(hgb_variants,pg.df,nomatch = NA) %>%  
  filter(PP>0.10) %>% 
  arrange(desc(PP)) %>% 
  dplyr::select(var,PP,FC,gene,qvalue,cor) %>%
  distinct()

mcv_variants <- as.data.table(mcv_variants)
setkey(setDT(mcv_variants), seqnames, start, end)

mcv_atac_cor<- foverlaps(mcv_variants,pg.df) %>%  
  filter(PP>0.10) %>% 
  arrange(desc(PP)) %>% 
  dplyr::select(var,PP,FC,gene,qvalue,cor) %>%
  distinct()

# MotifbreakR overlap
# HGB
motifbreakr <- readRDS("../../../singlecell_bloodtraits/data/motifbreakR/alltraits.mbreaker.withPPs.rds")

hgb_atac_cor$tomerge <- paste0("chr",gsub("_",":",hgb_atac_cor$var))
hgb_atac_cor_motifs <- left_join(hgb_atac_cor,motifbreakr[motifbreakr$trait %in% "HGB",],
                                 by=c("tomerge"="SNP")) %>% dplyr::select(c(1:6),geneSymbol) %>% dplyr::rename(motif=geneSymbol)

# MCV
mcv_atac_cor$tomerge <- paste0("chr",gsub("_",":",mcv_atac_cor$var))
mcv_atac_cor_motifs <- left_join(mcv_atac_cor,motifbreakr[motifbreakr$trait %in% "MCV",],
                             by=c("tomerge"="SNP")) %>% dplyr::select(c(1:6),geneSymbol) %>% 
  dplyr::rename(motif=geneSymbol) %>% distinct()


# PChic HGB variants
pchic <- readRDS("/Users/erikbao/Documents/GitHub/singlecell_bloodtraits/data/pchic/pchic.rds")
pchic.gr <- GRanges(pchic)
grch38.pc <- grch38 %>%
  dplyr::filter(biotype == "protein_coding")
pchic.gr <- pchic.gr[pchic.gr$Gene %in% grch38.pc$symbol,]

erytraits <- c("HCT", "HGB", "MCH", "MCHC", "MCV", "RBC_COUNT")
pchic.gr.RBC <- pchic.gr[pchic.gr$CellType == "Ery" & pchic.gr$Trait %in% erytraits,]
pchic.df.RBC <- as.data.frame(pchic.gr.RBC)
pchic.df.RBC$tomerge <- paste(pchic.df.RBC$seqnames,pchic.df.RBC$variantPos,sep=":")

hgb_atac_cor_motifs$tomerge <- str_split_fixed(paste0("chr",hgb_atac_cor_motifs$var),"_",3)[,1]
hgb_atac_cor_motifs_pchic <- left_join(hgb_atac_cor_motifs, pchic.df.RBC[pchic.df.RBC$Trait %in% "HGB",],by=c("tomerge"))  %>% 
  dplyr::select(c(1:7),Gene) %>% dplyr::rename(pchic_gene = Gene) %>% distinct()

mcv_atac_cor_motifs$tomerge <- str_split_fixed(paste0("chr",mcv_atac_cor_motifs$var),"_",3)[,1]
mcv_atac_cor_motifs_pchic <- left_join(mcv_atac_cor_motifs, pchic.df.RBC[pchic.df.RBC$Trait %in% "MCV",],by=c("tomerge"))  %>% 
  dplyr::select(c(1:7),Gene) %>% dplyr::rename(pchic_gene = Gene) %>% distinct()

write.table(hgb_atac_cor_motifs_pchic,"../../data/hgb_hct_variants_goodpeaks.tsv",sep = "\t", quote = FALSE, col.names = T, row.names = FALSE)
write.table(mcv_atac_cor_motifs_pchic,"../../data/mcv_rbc_mchc_mch_variants_goodpeaks.tsv",sep = "\t", quote = FALSE, col.names = T, row.names = FALSE)
