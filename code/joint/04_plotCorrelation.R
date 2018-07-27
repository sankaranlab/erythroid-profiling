library(data.table)
library(motifmatchr)
library(dplyr)
library(BuenColors)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(chromVARmotifs)
library(SummarizedExperiment)

allCor <- data.frame(fread(paste0("zcat < ", "../../processed/peakGeneCorrelations_1Mb_all.tsv.gz")))
allCor$FDR <- p.adjust(allCor$V6, method = "fdr")

plot_cor <-allCor %>% filter(FDR < 0.01)
ggplot(plot_cor, aes(x = V5)) + geom_density() +
  pretty_plot() + L_border()

plot_cor$activator <- plot_cor$V5 > 0
plot_cor$repressor <- plot_cor$V5 < 0

plot_cor %>% group_by(V1, V2, V3) %>%
  summarise(mA = mean(activator), 
            mR = mean(repressor)) %>% 
  filter(mA > 0 & mR > 0) %>% dim() -> n_both

plot_cor %>% group_by(V1, V2, V3) %>%
  summarise(mA = mean(activator), 
            mR = mean(repressor)) %>% 
  filter(mA > 0 ) %>% dim() -> n_activator

plot_cor %>% group_by(V1, V2, V3) %>%
  summarise(mA = mean(activator), 
            mR = mean(repressor)) %>% 
  filter(mR > 0 ) %>% dim() -> n_repressor


# Import peak set
peaksdf <- data.frame(fread("../../data/corces/panHeme.bed"))
peaks_gr <- makeGRangesFromDataFrame(peaksdf, seqnames.field = "V1", start.field = "V2", end.field = "V3")

propA <- n_activator[1]/dim(peaksdf)[1]
propR <- n_repressor[1]/dim(peaksdf)[1]

propA*propR*dim(peaksdf)[1]
n_both[1]

# If the expected < observed, certain enhancers more critical

# Do a motif enrichment here via fisher tests
data("human_pwms_v2")
matchm <- motifmatchr::matchMotifs(human_pwms_v2, peaks_gr, genome = BSgenome.Hsapiens.UCSC.hg19)

activate_gr <- plot_cor %>% filter(activator & !repressor) %>% distinct() %>% 
  makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")
a_ov <- findOverlaps(activate_gr, peaks_gr)

repress_gr <- plot_cor %>% filter(!activator & repressor) %>% distinct() %>% 
  makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")
r_ov <- findOverlaps(repress_gr, peaks_gr)

bdf <- data.frame(data.matrix(assays(matchm)[["motifMatches"]]))
bdf$activate <- 1:dim(bdf)[1] %in% subjectHits(a_ov)
bdf$repress <- 1:dim(bdf)[1] %in% subjectHits(r_ov)

enrichOut <- lapply(names(human_pwms_v2), function(tf){
  OR_P <- function(type){
    tt <- bdf[,c(tf, type)]
    ft <- fisher.test(matrix(c( sum(tt[,1] & tt[,2]), sum(!tt[,1] & tt[,2]), sum(tt[,1] & !tt[,2]), sum(!tt[,1] & !tt[,2])),ncol =2))
    return(c(ft$estimate, ft$p.value))
  }
  a <- OR_P("activate")
  r <- OR_P("repress")
  data.frame(ORactivate = a[1], ORprepress = r[1], tf)
}) %>% rbindlist() %>% data.frame()


colnames(enrichOut) <- c("OR_activate", "OR_repress", "TF")

p1 <- ggplot(enrichOut, aes(x = OR_activate, y = OR_repress, group = TF)) + geom_point()
plotly::ggplotly(p1)

