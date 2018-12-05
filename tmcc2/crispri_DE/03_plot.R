library(BuenColors)
library(data.table)
library(dplyr)

neighbors <- c("NFASC", "CNTN2", "RBBP5", "DSTYK", "NUAK2", "KLHDC8A", "LEMD1", "CDK18", "TMEM81")

df <- fread("output/RNAseq-g1_NC_Padj01.tsv") %>% data.frame()
df$color <- ifelse(df$gene == "TMCC2", "TMCC2", ifelse(df$gene %in% neighbors, "Neighbor", "other"))

ggplot(df %>% filter(baseMean > 20) %>% arrange(desc(color)), aes(x  =-1* log2FoldChange, y = -log10(pvalue), color = color)) +
  geom_point() +
  labs(x = "log2FC CRISPRi/Non-targeting", y = "-log10 P-value") + pretty_plot(fontsize = 8) + L_border() + 
  scale_color_manual(values = c("TMCC2" = "firebrick", "Neighbor" = "dodgerblue2", "other" = "black"))
  theme(legend.position = "none") 



