library(BuenColors)
library(dplyr)
library(data.table)

lapply(paste0("P", as.character(1:8)), function(population){
  data.frame(population, read.table(paste0("../kernel/", population, ".roc.out"))[,c(3,4)])
}) %>% rbindlist() %>% data.frame() -> df

colnames(df) <- c("pop", "AUROC", "AUPRC")
mdf <- reshape2::melt(df)

p1 <- ggplot(mdf, aes(x = pop, y = value, fill= variable)) +
  geom_bar(position=position_dodge(), stat="identity",
           color = "black", width = 0.6) +
  scale_fill_manual(values = c("dodgerblue3", "firebrick")) + 
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "Population", y = "Statistic") +
  scale_y_continuous( expand = c(0.1, 0)) +
  theme(legend.position = "bottom")

cowplot::ggsave(p1, file = "../../plots/ROC-PR.pdf", width = 3, height= 2.5)