library(dplyr)
library(diffloop)
library(GenomicRanges)
library(data.table)
library(BuenColors)

"%ni%" <- Negate("%in%")

x01_coding_gr <- bedToGRanges("../../data/annotations/Coding_UCSC.bed")
x02_promoter_gr <- bedToGRanges("../../data/annotations/Promoter_UCSC.fixed.bed")
x03_utr_gr <- bedToGRanges("../../data/annotations/UTR_3_UCSC.bed")
x04_SE_gr <- bedToGRanges("../../data/annotations/SE.bed")
x05_intron_gr <- bedToGRanges("../../data/annotations/Intron_UCSC.bed")

# Import peaks 
gr_t <- bedToGRanges("../../data/ATAC_data/ery_only.bed")

# Do overlaps
ov_1 <- findOverlaps(gr_t, x01_coding_gr)
ov_2 <- findOverlaps(gr_t, x02_promoter_gr)
ov_3 <- findOverlaps(gr_t, x03_utr_gr)
ov_4 <- findOverlaps(gr_t, x04_SE_gr)
ov_5 <- findOverlaps(gr_t, x05_intron_gr)

# Classify each accessibility peak
classAll <- ifelse(1:length(gr_t) %in% queryHits(ov_4), "SE",
                  ifelse(1:length(gr_t) %in% queryHits(ov_2), "promoter",
                         ifelse(1:length(gr_t) %in% queryHits(ov_3), "utr",
                                ifelse(1:length(gr_t) %in% queryHits(ov_5), "intron",
                                       ifelse(1:length(gr_t) %in% queryHits(ov_1), "coding", "distal")))))

# Import differentially accessible peaks
lapply(list.files("../../downloads/ATAC_DESeq2/", full.names = TRUE), function(x){
  data.frame(fread(x, header = TRUE, stringsAsFactors = FALSE))[,c(1,2,3)]
}) %>% rbindlist() %>% data.frame() %>% distinct() -> diffPeakDf


gr_t <- makeGRangesFromDataFrame(diffPeakDf)

# Do overlaps
ov_1 <- findOverlaps(gr_t, x01_coding_gr)
ov_2 <- findOverlaps(gr_t, x02_promoter_gr)
ov_3 <- findOverlaps(gr_t, x03_utr_gr)
ov_4 <- findOverlaps(gr_t, x04_SE_gr)
ov_5 <- findOverlaps(gr_t, x05_intron_gr)

classD <- ifelse(1:length(gr_t) %in% queryHits(ov_4), "SE",
                  ifelse(1:length(gr_t) %in% queryHits(ov_2), "promoter",
                         ifelse(1:length(gr_t) %in% queryHits(ov_3), "utr",
                                ifelse(1:length(gr_t) %in% queryHits(ov_5), "intron",
                                       ifelse(1:length(gr_t) %in% queryHits(ov_1), "coding", "distal")))))


data.frame(all=table(classAll)*100/sum(table(classAll)),
           differential=table(classD)*100/sum(table(classD)))

