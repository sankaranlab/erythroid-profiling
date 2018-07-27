library(BuenColors)
library(dplyr)
library(data.table)

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")
counts <- data.matrix(data.frame(fread("../../data/corces/atac_combinedPopulations/panHeme.counts.tsv")))

makeFactorPlot <- function(factor){
  # Write a function to import and normalize the Tn5 density
  getInsertScores <- function(population, factor){
    #idxstats <- data.frame(fread(paste0("../../data/ATAC_data/idxstats/", population, ".bam.idxstats")))
    #total <- sum(idxstats$V3)
    total <- sum(counts[,population])
    tn5 <- readRDS(paste0("../../data/ATAC_data/tn5_density/", population, "-", factor, ".tag.rds"))
    outdf <- data.frame(
      population = population, 
      position = tn5$pos,
      value = tn5$value * (100000/total)
    )
    
  }
  
  
  gata1inserts <- rbind(
    getInsertScores("P1", factor),
    getInsertScores("P5", factor),
    getInsertScores("P8", factor))
  
  p1 <- ggplot(gata1inserts, aes(x = position, y = value, color = population, group = population)) +
    geom_line(size = 0.25) + scale_color_manual(values = eryth_color_maps) + 
    pretty_plot(fontsize = 8) + L_border() + labs(x = paste0("Position relative to center of ",factor," Motif"),
                                                  y = "Normalized Tn5 Insertions", color = "") +
    theme(legend.position = "none") 
  factor
  cowplot::ggsave(p1, file = paste0("../../plots/",factor,"-inserts.pdf"), width = 2, height = 2)

}

makeFactorPlot("CTCF")
makeFactorPlot("KLF14")
makeFactorPlot("GATA1")
makeFactorPlot("IRF1")
makeFactorPlot("FOSL1")
makeFactorPlot("NFE2")
makeFactorPlot("SPIB")
