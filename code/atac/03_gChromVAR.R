library(chromVAR)
library(gchromVAR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(tidyverse)
library(SummarizedExperiment)
library(Matrix)
library(BuenColors)
library(chromVARmotifs)
library(cowplot)
library(diffloop)
library(preprocessCore)
library(stringr)

"%ni%" <- Negate("%in%")


# Import and run run run 
peaksdf <- fread("../../data/corces/panHeme.bed")
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(data.frame(fread("../../data/corces/atac_combinedPopulations/panHeme.counts.tsv")))

# Create objects for g-chromVAR
SE <- SummarizedExperiment(assays = list(counts = counts),
                           rowData = peaks, 
                           colData = DataFrame(names = colnames(counts)))
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
ukbb <- importBedScore(rowRanges(SE), list.files("../../data/FMsnps/", full.names = TRUE, pattern = "*.bed$"))

# Run g-chromVAR
bg <- getBackgroundPeaks(SE)
dev <- computeWeightedDeviations(SE, ukbb, background_peaks = bg)
outdf <- assays(dev)[["z"]]
rownames(outdf) <- gsub("_PP001_betas", "", rownames(outdf))
mdf <- reshape2::melt(outdf)
mdf$Var2 <- factor(as.character(mdf$Var2), levels = c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8))))

ggplot(mdf[mdf$Var1 %in% c( "HGB", "MCV", "MCH", "MCHC", "RETIC_COUNT",
                            "MEAN_RETIC_VOL", "RBC_COUNT","PLT_COUNT") , ],
       aes(x = Var2, y = value, color = Var1, group = Var1)) +
  geom_point() + geom_line() + pretty_plot(fontsize = 8) + L_border() + 
  labs(x = "", y = "g-chromVAR Zscore", color = "") + theme(legend.position = "right") 

ggsave(p1, filename = "../../plots/gchromVAR_selected.pdf", width = 3.5, height = 1.5)


# g-chromVAR using P1-P8 only --------------------------------------------------------------

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
erytraits <- c("HCT", "HGB", "MCH", "MCHC", "MCV", "MEAN_RETIC_VOL", "RBC_COUNT", "RETIC_COUNT")

# Run g-chromVAR
bg <- getBackgroundPeaks(SE)
dev <- computeWeightedDeviations(SE, ukbb, background_peaks = bg)
outdf <- assays(dev)[["z"]]
rownames(outdf) <- gsub("_PP001_betas", "", rownames(outdf))
mdf <- reshape2::melt(outdf)
mdf$Var2 <- factor(as.character(mdf$Var2), levels = c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8))))

traits_to_plot <- c("HGB","HCT","MCV","MCHC","MCH","RBC_COUNT") # Traits endorsed by Leif Ludwig
p2 <-ggplot(mdf[mdf$Var1 %in% traits_to_plot , ],
             aes(x = Var2, y = value, color = Var1, group = Var1)) +
  geom_point(size=2) + geom_line(size=1) + 
  pretty_plot(fontsize = 10) + L_border() +
  scale_color_manual(values=c(jdb_palette("Rushmore"),jdb_palette("flame_light")[c(5)])) +
  labs(x = "", y = "g-chromVAR Zscore", color = "") + theme(legend.position = "bottom") 

ggsave(p2, filename = "../../plots/gchromVAR_eryonly_MCV_MCH_MCHC_RBC.pdf", width = 3.5, height = 3.5)
