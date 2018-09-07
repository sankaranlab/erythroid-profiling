library(dplyr)
library(BuenColors)
library(data.table)
library(ggbeeswarm)

lapply(list.files("../../data/ATAC_data/summits/"), function(x){
  df <- fread(paste0("../../data/ATAC_data/summits/", x), header = TRUE)
  y <- stringr::str_split_fixed(x, "_", 4)
  data.frame(Sample = gsub("_ATAC_summits.bed", "", paste0(y[1,3], "-", y[1,4])),
    Population = stringr::str_split_fixed(x, "_", 4)[1,3], nPeaks = dim(df)[1])
}) %>% rbindlist() %>% data.frame(stringsAsFactors = FALSE) %>% arrange(Sample) -> peakDF

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")

p1  <- ggplot(peakDF, aes(x=Population, y=(nPeaks), color = "X", fill = Population)) + pretty_plot(fontsize = 8) +
  stat_summary(fun.y="mean", geom="bar") + scale_color_manual(values = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4) +  geom_quasirandom(size = 0.5) +L_border() +
  scale_fill_manual(values = eryth_color_maps) + labs(x = "", y = "# Accessibility Peaks") +
  theme(legend.position = "none")

cowplot::ggsave(p1, file = "../../plots/nPeaks_ATAC.pdf", width = 3, height = 2)
