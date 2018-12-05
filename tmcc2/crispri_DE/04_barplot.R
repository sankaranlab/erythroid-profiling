library(dplyr)
library(BuenColors)

# Import raw expression values
raw <- read.table("output/TMCC2-RNAseq_rawGeneCounts.tsv", header = TRUE)
RNA.counts <- raw[,1:7]
cpm <- round(sweep(RNA.counts, 2, colSums(RNA.counts), FUN="/") * 1000000, 1)
rownames(cpm) <- raw[,8]

cpm[intersect(c("NFASC", "CNTN2", "RBBP5", "DSTYK", "NUAK2", "KLHDC8A", "LEMD1", "CDK18", "TMEM81"), rownames(cpm)),]

# Select for the non-zero
pick <- c("RBBP5", "DSTYK", "KLHDC8A")


data.frame(t(cpm["TMCC2",])) %>% 
  mutate( V1 = c("Control", "Control", "Control", "Enhancer", "Enhancer", "Enhancer", "Enhancer"),
          V2 = c("Control", "Control", "Control", "Guide1", "Guide1", "Guide2", "Guide2")) %>%
  group_by(V1, V2) %>% summarise(mean = mean(TMCC2),
                                       sd = sd(TMCC2)) %>% data.frame() -> tmcc2

cpm[pick,] %>% t() %>% data.frame() %>%
  mutate( V1 = c("Control", "Control", "Control", "Enhancer", "Enhancer", "Enhancer", "Enhancer")) %>%
  reshape2::melt(id.vars = "V1") %>% 
  group_by(V1, variable) %>% summarise(mean = mean(value),
                                       sd = sd(value)) %>% data.frame() -> neighbor_df


p0 <- ggplot(tmcc2, aes(x = V2, y = mean, fill = V1)) +
  pretty_plot(fontsize = 8) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge(), alpha = 0.3) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("Control" = "black", "Enhancer" = "firebrick")) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  #  geom_beeswarm(inherit.aes = FALSE, data = points, aes(x = V1, y = V3)) +
  labs(x = "", y = "TMCC2 expression (cpm)", fill = "gRNA") + theme(legend.position = c(0.8, 0.8))


p1 <- ggplot(neighbor_df, aes(x = variable, y = mean, fill = V1)) +
  pretty_plot(fontsize = 8) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge(), alpha = 0.3) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("Control" = "black", "Enhancer" = "firebrick")) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  #  geom_beeswarm(inherit.aes = FALSE, data = points, aes(x = V1, y = V3)) +
  labs(x = "", y = "RNA-seq expression (cpm)", fill = "gRNA") + theme(legend.position = c(0.8, 0.8))

# Do the volcano
neighbors <- c("NFASC", "CNTN2", "RBBP5", "DSTYK", "NUAK2", "KLHDC8A", "LEMD1", "CDK18", "TMEM81")

df <- fread("output/RNAseq-g-pooled_NC_Padj01.tsv") %>% data.frame()
df$color <- ifelse(df$gene == "TMCC2", "TMCC2", ifelse(df$gene %in% neighbors, "Neighbor", "other"))

pV <- ggplot(df %>% filter(baseMean > 20) %>% arrange(desc(color)), aes(x  =-1* log2FoldChange, y = -log10(pvalue), color = color)) +
  geom_point(size = 0.5) +
  labs(x = "log2FC CRISPRi/Non-targeting", y = "-log10 P-value") + pretty_plot(fontsize = 8) + L_border() + 
  scale_color_manual(values = c("TMCC2" = "firebrick", "Neighbor" = "dodgerblue2", "other" = "black")) +
  theme(legend.position = "none") 


cowplot::ggsave(cowplot::plot_grid(p0, p1, pV, nrow = 1), 
                file = "three_summary.pdf", width = 6.3, height = 2.1)


