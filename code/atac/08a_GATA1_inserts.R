library(BuenColors)
library(dplyr)
library(data.table)

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")

# Write a function to import and normalize the Tn5 density
getInsertScores <- function(population, factor){
  idxstats <- data.frame(fread(paste0("../../data/ATAC_data/idxstats/", population, ".bam.idxstats")))
  total <- sum(idxstats$V3)
  tn5 <- readRDS(paste0("../../data/ATAC_data/tn5_density/", population, "-", factor, ".tag.rds"))
  outdf <- data.frame(
    population = population, 
    position = tn5$pos,
    value = tn5$value * (100000/total)
  )
  
}

gata1inserts <- rbind(
  getInsertScores("P1", "GATA1"),
  getInsertScores("P5", "GATA1"),
  getInsertScores("P8", "GATA1"))

p1 <- ggplot(gata1inserts, aes(x = position, y = value, color = population, group = population)) +
  geom_line(size = 0.25) + scale_color_manual(values = eryth_color_maps) + 
  pretty_plot(fontsize = 8) + L_border() + labs(x = "Position relative to center of GATA1 Motif",
                                    y = "Normalized Tn5 Insertions", color = "") +
  theme(legend.position = "none") 

cowplot::ggsave(p1, file = "../../plots/GATA1-inserts.pdf", width = 2, height = 2)