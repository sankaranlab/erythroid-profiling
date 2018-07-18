library(data.table)
library(motifmatchr)
library(dplyr)
library(BuenColors)

allCor <- data.frame(fread(paste0("zcat < ", "../../processed/peakGeneCorrelations_1Mb_all.tsv.gz")))
allCor$FDR <- p.adjust(allCor$V6, method = "fdr")

plot_cor <-allCor %>% filter(FDR < 0.05)
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

propA <- n_activator[1]/dim(peaksdf)[1]
propR <- n_repressor[1]/dim(peaksdf)[1]

propA*propR*dim(peaksdf)[1]
n_both[1]

# If the expected < observed, certain enhancers more critical


# Do a motif enrichment here via fisher tests


