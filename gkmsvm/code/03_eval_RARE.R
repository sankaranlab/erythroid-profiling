library(BuenColors)
library(dplyr)

importOne <- function(pop) read.table(paste0("../variantPredictions/",pop,"_mendelianRare.out"))[,2]

pred_df <- data.frame(sapply(paste0("P", as.character(1:8)), importOne))
pred_df$variant <- read.table(paste0("../variantPredictions/P1_mendelianRare.out"))[,1]

mdf <- reshape2::melt(pred_df, id.var = "variant")

mdf %>% 
  group_by(variant) %>% mutate(minMaxColor = (value - min(value))/(max(value) - min(value))) %>%
  arrange(variant, desc(value)) -> hemeGrid

p1 <- ggplot(hemeGrid, aes(x = variant, y = variable, fill = value)) +
  geom_tile( color = "black") +
  scale_fill_gradientn(colors = jdb_palette("brewer_spectra")) +
  labs(x = "", y = "", fill = "Enrichment") + pretty_plot(fontsize = 8) +
  theme(legend.position = "bottom") + L_border() + 
  theme(axis.text.x=element_text(angle=45, hjust=1))

