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
"%ni%" <- Negate("%in%")

# Identifying g-chromVAR trends -------------------------------------------

# Import and run run run 
peaksdf <- fread("../../data/ATAC_data/ery_only.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")

ATAC.counts <- data.matrix(data.frame(data.table::fread("../../data/ATAC_data/ery_only.counts.tsv")))
meta <- stringr::str_split_fixed(colnames(ATAC.counts), "_", 4)
pops.ordered <- c( paste0("P", as.character(1:8)))

sapply(pops.ordered, function(pop){
  rowSums(ATAC.counts[,meta[,3] == pop])
}) -> counts

# Create objects for g-chromVAR
SE <- SummarizedExperiment(assays = list(counts = counts),
                           rowData = peaks, 
                           colData = DataFrame(names = colnames(counts)))
SE <- filterPeaks(SE)
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
ukbb <- importBedScore(rowRanges(SE), list.files("../../data/FMsnps/", full.names = TRUE, pattern = "*.bed$"))

# Tangent-- find peaks with high PP + accessibility 
cpm <- sweep(counts, 2, colSums(counts), FUN="/") * 1000000

# Take only peaks with a cumulative PP of 0.1 across all erythroid traits
erytraits <- c("HCT", "HGB", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL", "RBC_COUNT", "RETIC_COUNT")
keepPeaks <- rowSums(assays(ukbb)[["weights"]][,paste0(erytraits, "_PP001_betas")]) > 0.5
mega_df <- data.frame(
  rowRanges(SE)[keepPeaks],
  data.matrix(cpm[keepPeaks,]), # fix here with CPM
  data.matrix(assays(ukbb)[["weights"]][keepPeaks,paste0(erytraits, "_PP001_betas")]))
colnames(mega_df)[15:ncol(mega_df)] <- erytraits

mega_df$meanP2_4 <- apply(mega_df[,paste0("P",seq(2,4))],1,mean)
mega_df$meanP5_6 <- apply(mega_df[,paste0("P",seq(5,6))],1,mean)
mega_df$FC <- log2(mega_df$meanP5_6  / mega_df$meanP2_4)

# Overlap with finemapped variants
# CS.gr <- readRDS("../../data/UKBB_BC_v3_VEPannotations.rds")
CS.df <- CS.gr %>% as.data.frame() %>% as.data.table
setkey(setDT(CS.df), seqnames, start, end)

mega_df<-as.data.table(mega_df)
names(mega_df) <- c("chrom","j.start","j.end",names(mega_df)[4:ncol(mega_df)])
setkey(setDT(mega_df), chrom, j.start, j.end)

traits_to_plot <- c("HGB","HCT")
hgb_variants <- foverlaps(CS.df,mega_df, nomatch = 0)  %>%
  dplyr::filter(trait %in% traits_to_plot) %>% arrange(desc(PP)) %>% 
  dplyr::select(var,trait, seqnames, start,end,PP,meanP2_4,meanP5_6,FC) %>%
  distinct(var,.keep_all=TRUE)

# HGB/HCT Volcano plot
p1 <- ggplot(hgb_variants, aes(x=FC, y=PP)) +
  geom_point() +
  pretty_plot() + 
  L_border() + 
  geom_vline(xintercept = 0, linetype = 2)
ggExtra::ggMarginal(p1, type = "boxplot",margins="x")

plot(density(all_variants$FC))

# MCV/MCH Volcano plot
traits_to_plot <- c("MCV","MCH")
mcv_variants <- foverlaps(CS.df,mega_df, nomatch = 0)  %>%
  dplyr::filter(trait %in% traits_to_plot) %>% arrange(desc(PP)) %>% 
  dplyr::select(var,trait, seqnames, start,end,PP,meanP2_4,meanP5_6,FC) %>%
  distinct(var,.keep_all=TRUE)

# Volcano plot
ggplot(mcv_variants, aes(x=FC, y=PP)) +
  geom_point() +
  pretty_plot() + 
  L_border() + 
  geom_vline(xintercept = 0, linetype = 2)

plot(density(all_variants$FC))

exdf <- read.table("../../../singlecell_bloodtraits/figures/revisions/exclude_list_revised.txt", header = FALSE, stringsAsFactors = FALSE)[,1]

# Top HGB peaks
mega_df %>% filter(meanP2_4>meanP5_6) %>% filter(HGB>0.10) %>% arrange(desc(HGB)) -> HGB_peaks

pg.df <- fread(paste0("zcat < ", "/Users/erikbao/Documents/GitHub/singlecell_bloodtraits/data/bulk/peakGeneCorrelation.tsv.gz"))
names(pg.df) <- c("chrom","j.start","j.end","gene","cor","pvalue")
pg.df$qvalue <- qvalue(pg.df$pvalue)$qvalues
pg.df <- pg.df %>%
  dplyr::filter(qvalue < 0.001)

AR_cor <- data.frame()
HGB_peaks <- as.data.frame(HGB_peaks)
for (i in 1:nrow(HGB_peaks)){
  AR_cor <- rbind(AR_cor,pg.df[pg.df$chrom %in% HGB_peaks[i,"seqnames"] & 
                                 pg.df$j.start <= HGB_peaks[i,"start"]  &
                                 pg.df$j.end >= HGB_peaks[i,"end"] ,])
}

# Overlap with finemapped variants
CS.gr <- readRDS("../../data/UKBB_BC_v3_VEPannotations.rds")
CS.df <- CS.gr %>% as.data.frame() %>% as.data.table
setkey(setDT(CS.df), seqnames, start, end)

mega_df %>% filter(meanP2_4>meanP5_6) %>% filter(HGB>0.10) %>% arrange(desc(HGB)) -> HGB_peaks
HGB_peaks<-as.data.table(HGB_peaks)
names(HGB_peaks) <- c("chrom","j.start","j.end",names(HGB_peaks)[4:ncol(HGB_peaks)])
setkey(setDT(HGB_peaks), chrom, j.start, j.end)

HGB_variants <- foverlaps(CS.df,HGB_peaks, nomatch = 0)  %>%
  dplyr::filter(PP > 0.1 & trait == "HGB") %>% arrange(desc(PP)) %>% dplyr::select(var,seqnames, start,end,PP,meanP2_4,meanP5_6,FC)

setkey(setDT(HGB_variants), seqnames, start, end)
setkey(setDT(pg.df), chrom, j.start, j.end)
hgb_atac_cor<- foverlaps(HGB_variants,pg.df) %>%  arrange(desc(PP)) %>% 
  dplyr::select(var,seqnames, start,end,PP,meanP2_4,meanP5_6,FC,gene,qvalue)

# PChic HGB variants
pchic <- readRDS("/Users/erikbao/Documents/GitHub/singlecell_bloodtraits/data/pchic/pchic.rds")
pchic.gr <- GRanges(pchic)
grch38.pc <- grch38 %>%
  dplyr::filter(biotype == "protein_coding")
pchic.gr <- pchic.gr[pchic.gr$Gene %in% grch38.pc$symbol,]

erytraits <- c("HCT", "HGB", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL", "RBC_COUNT", "RETIC_COUNT")
pchic.gr.RBC <- pchic.gr[pchic.gr$CellType == "Ery" & pchic.gr$PP > 0.1 & pchic.gr$Trait %in% "HGB",]
pchic.df.RBC <- as.data.frame(pchic.gr.RBC)
pchic.df.RBC$tomerge <- paste(pchic.df.RBC$seqnames,pchic.df.RBC$variantPos,sep=":")
hgb_atac_cor$tomerge <- paste(hgb_atac_cor$seqnames,hgb_atac_cor$start,sep=":")

pchic_atac_merged <- left_join(hgb_atac_cor, pchic.df.RBC,by=c("tomerge")) 

# MotifbreakR overlap
motifbreakr <- readRDS("../../../singlecell_bloodtraits/data/motifbreakR/alltraits.mbreaker.withPPs.rds")

hgb_atac_cor$tomerge <- paste0("chr",gsub("_",":",hgb_atac_cor$var))

hgb_atac_cor_motifs <- left_join(hgb_atac_cor,motifbreakr[motifbreakr$trait %in% "HGB",],
                                 by=c("tomerge"="SNP")) %>% dplyr::select(c(1:9),geneSymbol) %>% 
  rename(motif=geneSymbol)

write.table(hgb_atac_cor_motifs,"../../data/hgb_variants.tsv",sep = "\t", quote = FALSE, col.names = T, row.names = FALSE)

vars_reformatted <- paste0("chr",gsub("_",":",hgb_atac_cor$var))
hgb_motifs <- motifbreakr[motifbreakr$SNP %in% vars_reformatted,]
hgb_motifs %>% filter(PP>0.10,trait=="HGB",effect=="strong") 

# Top MCV peaks
mega_df %>% filter(meanP5_6>meanP2_4) %>% filter(MCV>0.10) %>%  arrange(desc(FC)) -> MCV_peaks

AR_cor <- data.frame()
for (i in 1:nrow(MCV_peaks)){
  AR_cor <- rbind(AR_cor,pg.df[pg.df$chrom %in% MCV_peaks[i,"seqnames"] & 
                                 pg.df$j.start <= MCV_peaks[i,"start"]  &
                                 pg.df$j.end >= MCV_peaks[i,"end"] ,])
}

MCV_peaks<-as.data.table(MCV_peaks)
setkey(setDT(MCV_peaks), seqnames, start, end)
MCV_variants <- foverlaps(CS.df,MCV_peaks, nomatch = 0)  %>%
  dplyr::filter(PP > 0.1 & trait == "MCV") %>% arrange(desc(PP)) %>% dplyr::select(var,seqnames, start,end,PP,meanP2_4,meanP5_6,FC)

setkey(setDT(MCV_variants), seqnames, start, end)
setkey(setDT(pg.df), chrom, j.start, j.end)
atac_cor<- foverlaps(MCV_variants,pg.df) %>%  arrange(desc(PP)) %>% 
  dplyr::select(var,seqnames, start,end,PP,meanP2_4,meanP5_6,FC,gene,qvalue)

# MotifbreakR overlap
atac_cor$tomerge <- paste0("chr",gsub("_",":",atac_cor$var))
atac_cor_motifs <- left_join(atac_cor,motifbreakr[motifbreakr$trait %in% "MCV",],
                                 by=c("tomerge"="SNP")) %>% dplyr::select(c(1:9),geneSymbol) %>% 
  rename(motif=geneSymbol) %>% distinct()

write.table(atac_cor_motifs,"../../data/mcv_variants.tsv",sep = "\t", quote = FALSE, col.names = T, row.names = FALSE)


vars_reformatted <- paste0("chr",gsub("_",":",MCV_variants$var))
mcv_motifs <- motifbreakr[motifbreakr$SNP %in% vars_reformatted,]
mcv_motifs %>% filter(PP>0.10,trait=="MCV",effect=="strong") %>% dplyr::select(geneSymbol)


# Read in PCHIC 
library(annotables)
pchic <- readRDS("/Users/erikbao/Documents/GitHub/singlecell_bloodtraits/data/pchic/pchic.rds")
pchic.gr <- GRanges(pchic)
grch38.pc <- grch38 %>%
  dplyr::filter(biotype == "protein_coding")
pchic.gr <- pchic.gr[pchic.gr$Gene %in% grch38.pc$symbol,] 
pchic.df <- pchic.gr %>% as.data.frame

pchic_overlap <- data.frame()
for (i in 1:nrow(HGB_peaks)){
  pchic_overlap <- rbind(AR_cor,pg.df[pchic.df$seqnames %in% HGB_peaks[i,"seqnames"] & 
                                        pchic.df$variantPos >= HGB_peaks[i,"start"]  &
                                        pchic.df$end >=HGB_peaks[i,"end"] ,])
}

# Top ATAC Peaks in P2-P4 ------------------------------------------------------
# Import peaks
bed <- data.frame(fread("../../data/ATAC_data/ery_only.bed", header = FALSE))
colnames(bed) <- c("chr", "start", "end")
peaks_gr <- makeGRangesFromDataFrame(bed)

# Counts
ATAC.counts <- data.matrix(data.frame(data.table::fread("../../data/ATAC_data/ery_only.counts.tsv")))
meta <- stringr::str_split_fixed(colnames(ATAC.counts), "_", 4)
pops.ordered <- c( paste0("P", as.character(1:8)))

# Filter for only peaks where it's at least in top nth percentile for P2, P3, and/or P4
sapply(pops.ordered, function(pop){
  rowSums(ATAC.counts[,meta[,3] == pop])
}) -> counts.df.raw

counts.df <- counts.df.raw
n = 0.7
counts.df[1:dim(counts.df.raw)[1],1:dim(counts.df.raw)[2]] <- normalize.quantiles(as.matrix(counts.df.raw))
keep <- apply(counts.df[,c("P2","P3","P4")],1,max) > mean(apply(counts.df,2,function(x) {quantile(x,n)}))
counts.df.raw <- counts.df.raw[keep,]

ATAC.cpm <- round(sweep(counts.df.raw, 2, colSums(counts.df.raw), FUN="/") * 1000000, 1)
ATAC.cpm.log2 <- log2(ATAC.cpm+1)

# Subset to only good peaks
peaks.gr <- peaks_gr[keep,]

# Overlap strong P2-P4 peaks with FM variants -----------------------------

#CS.gr <- readRDS("../../../singlecell_bloodtraits/data/Finemap/UKBB_BC_v3_VEPannotations.rds")

ov1 <- findOverlaps(CS.gr, peaks.gr)
P2_P4_df <- data.frame(CS.gr[queryHits(ov1)], round(ATAC.cpm[subjectHits(ov1),], 1))
P2_P4_df <- P2_P4_df[P2_P4_df$PP > 0.10,]
P2_P4_df[P2_P4_df$trait %in% c("HCT"),]

library(motifbreakR)
motifbreakr <- readRDS("../../../singlecell_bloodtraits/data/motifbreakR/alltraits.mbreaker.withPPs.rds")

# Motifs disrupted by HGB variants in P2-P4
P2_P4_df[P2_P4_df$trait %in% "HGB",]$var -> k6_vars

k6_vars_reformatted <- paste0("chr",gsub("_",":",k6_vars))
k6_hgb_motifs <- motifbreakr[motifbreakr$SNP %in% k6_vars_reformatted,]
k6_hgb_motifs %>% filter(PP>0.10,trait=="HGB",effect=="strong") 
k6_hgb_motifs %>% filter(PP>0.10,trait=="HGB",effect=="strong") -> output
write.table(output$geneSymbol,file="../../data/FUMA/HGB_analysis/k6_0.7_mbreaker.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Overlap strong P5-6 peaks with non-HGB/HCT variants --------------------------------

# Remove "weak" peaks that aren't in top n% for P2-P4 compared to all peaks
sapply(pops.ordered, function(pop){
  rowSums(ATAC.counts[,meta[,3] == pop])
}) -> counts.df.raw

counts.df <- counts.df.raw
n = 0.7
counts.df[1:dim(counts.df.raw)[1],1:dim(counts.df.raw)[2]] <- normalize.quantiles(as.matrix(counts.df.raw))
keep <- apply(counts.df[,c("P5","P6")],1,max) > mean(apply(counts.df,2,function(x) {quantile(x,n)}))
counts.df.raw <- counts.df.raw[keep,]

ATAC.cpm <- round(sweep(counts.df.raw, 2, colSums(counts.df.raw), FUN="/") * 1000000, 1)
ATAC.cpm.log2 <- log2(ATAC.cpm+1)

# Subset to only good peaks
P5_7_peaks.gr <- peaks_gr[keep,]

# Overlap
traits_of_interest <- c("MCV","MCHC")
ov1 <- findOverlaps(CS.gr, P5_7_peaks.gr)
P5_P6_df <- data.frame(CS.gr[queryHits(ov1)], round(ATAC.cpm[subjectHits(ov1),], 1))
P5_P6_df <- P5_P6_df[P5_P6_df$PP > 0.10,]
P5_P6_df[P5_P6_df$trait %in% traits_of_interest ,]

# Motifs disrupted by variants in P5-6
P5_P6_df[P5_P6_df$trait %in% traits_of_interest,]$var -> p56_vars

p56_vars_reformmatted <- paste0("chr",gsub("_",":",p56_vars))
k6_hgb_motifs <- motifbreakr[motifbreakr$SNP %in% p56_vars_reformmatted,]
k6_hgb_motifs %>% filter(PP>0.10,trait=="HGB",effect=="strong") 
k6_hgb_motifs %>% filter(PP>0.10,trait=="HGB",effect=="strong") -> output
write.table(output$geneSymbol,file="../../data/FUMA/HGB_analysis/k6_0.7_mbreaker.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)


# Look at HGB atac correlation peaks --------------------------------------

# Read output 
erycor.df <- fread(paste0("zcat < ","../../processed/peakGeneCorrelations_1Mb_all.tsv.gz"))
colnames(erycor.df) <- c("seqnames","start","end","gene","rho","p")
erycor.df$qvalue <- qvalue(erycor.df$p)$qvalues
erycor.df <- erycor.df %>%
  dplyr::filter(qvalue < 0.001)

# Calculate median number of unique peaks correlated with expression of each gene
temp <- erycor.df %>% group_by(gene) %>% summarize(n())
median(temp$`n()`)
plot(density(temp$`n()`))

# Calculate number of peaks associated with hemoglobin gene expression
# hgb_genes <- c("HBA1","HBA2","HBB","HBD","HBE1",
#                "HBG1","HBG2","HBM","HBQ1","HBZ")
hgb_genes <- c("HBA1","HBB")
temp2 <- erycor.df[erycor.df$gene %in% hgb_genes,] %>% group_by(gene) %>% summarize(n()) 
median(temp2$`n()`)


# chromVAR comparison P2-4,P5-6 ------------------------------------------------
dev <- readRDS(file="../../data/chromVAR/chromVAR_eryonly_dev.rds")
devscores <- read.table("../../data/chromVAR/chromVAR_eryonly_devscores.txt")
bagged <- chromVARxx::bagDeviations(dev, 0.8, "human")
bagged_devscores <- as.data.frame(deviationScores(bagged))

extract_top_subset <- function(df,subset=c("P2","P3","P4"),n=0.7){
  df[1:dim(df)[1],1:dim(df)[2]] <- normalize.quantiles(as.matrix(df))
  keep <- apply(df[,subset],1,max) > mean(apply(df,2,function(x) {quantile(x,n)}))
  df[keep,] %>% tibble::rownames_to_column() %>% as.tibble()
}

p2_p4_devscores <- extract_top_subset(devscores,subset=c("P2","P3","P4"),n=0.7) %>% arrange(desc(P4))
p2_p4_bagged_devscores <- extract_top_subset(df=bagged_devscores,subset=c("P2","P3","P4"),n=0.7) %>% arrange(desc(P4))

p56_devscores <- extract_top_subset(devscores,subset=c("P5","P6"),n=0.7) %>% arrange(desc(P5))
p56_bagged_devscores <- extract_top_subset(df=bagged_devscores,subset=c("P5","P6"),n=0.7) %>% arrange(desc(P5))

# Generate Heatmap of Hgb-correlated peaks across heme states  ---------------------

# Import peaks
bed <- data.frame(fread("../../data/corces/panHeme.bed", header = FALSE))
colnames(bed) <- c("chr", "start", "end")
peak_gr <- makeGRangesFromDataFrame(bed)

# Import raw expression values
raw <- data.frame(fread("../../data/corces/atac_combinedPopulations/panHeme.counts.tsv", header = TRUE))
ATAC.counts <- data.matrix(raw)
meta <- colnames(ATAC.counts)

# Normalize
cpm <- sweep(raw, 2, colSums(raw), FUN="/") * 1000000

hgb_peaks <- erycor.df[erycor.df$gene %in% hgb_genes,] %>% makeGRangesFromDataFrame
#hgb_peaks <- pg.df[pg.df$gene %in% hgb_genes,] %>% makeGRangesFromDataFrame

ov1 <- findOverlaps(hgb_peaks, peak_gr)
hgb_countsonly <- round(cpm[subjectHits(ov1),], 1)
hgb_countsonly <- cpm[sample(nrow(cpm),size=1000,replace=FALSE),]

ATAC.cpm.log2.all.mm.Hgb <- t(apply(hgb_countsonly, 1, function(x)(x-min(x))/(max(x)-min(x))))

# Heatmap(ATAC.cpm.log2.all.mm.Hgb, col=as.character(jdb_palette("solar_rojos",type="continuous")),
#         cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
#         row_names_gp = gpar(fontsize = 0),
#         column_names_gp = gpar(fontsize = 6),
#         split = km.cluster, show_heatmap_legend = FALSE,
#         name = "Peak\nAccessibility")

Heatmap(ATAC.cpm.log2.all.mm.Hgb,cluster_columns=F,show_row_names = F)
