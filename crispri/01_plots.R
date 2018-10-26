library(BuenColors)
library(ggbeeswarm)
library(dplyr)

points <- read.table("rawData-satish-oct25.txt")

points %>% group_by(V1) %>% summarise(mean = mean(V3),
                                          sd = sd(V3)) %>% data.frame() -> qPCR_df_plot
qPCR_df_plot$V4 <- c("Control", "enhancer", "enhancer", "promoter", "promoter")

p1 <- ggplot(qPCR_df_plot, aes(x = V1, y = mean, fill = V4)) +
  pretty_plot(fontsize = 8) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge(), alpha = 0.3) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("Control" = "black", "promoter" = "grey", "enhancer" = "firebrick")) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
#  geom_beeswarm(inherit.aes = FALSE, data = points, aes(x = V1, y = V3)) +
  labs(x = "", y = "TMCC2 Relative expression") + theme(legend.position = "none")


cowplot::ggsave(p1,
                file = "viewMe.pdf", width = 3.5, height = 1.9)
