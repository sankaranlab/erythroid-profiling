library(dplyr)
library(BuenColors)

# 1) qPCR data
qPCR_df <- data.frame(
  shRNA = c(rep("shluc", 3), rep("sh1", 3),  rep("sh2", 3)),
  val = c(1.29, 0.73, 1.06, 0.02, 0.15, 0.01, 0.03, 0.01, 0.01)
)

qPCR_df %>% group_by(shRNA) %>% summarise(mean = mean(val),
                                          sd = sd(val)) %>% data.frame() -> qPCR_df_plot
qPCR_df_plot$shRNA <- factor(as.character(qPCR_df_plot$shRNA), c("shluc", "sh1", "sh2"))

# make plot one
p1 <- ggplot(qPCR_df_plot, aes(x = shRNA, y = mean, fill = shRNA)) +
  pretty_plot(fontsize = 8) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("shluc" = "black", "sh2" = "orange3", "sh1" = "firebrick")) + L_border() +
  labs(x = "", y = "Relative Expression") + theme(legend.position = "none") + 
  scale_y_continuous(expand = c(0,0))

# 2) growth data
growth_df <- data.frame(
  shRNA = c(rep("shluc", 6), rep("sh1", 6),  rep("sh2", 6)),
  day = rep(c(rep("D10", 2), rep("D12", 2), rep("D15", 2)),3),
  val = c(1,1,1.47,1.42,2.23,2.26,1,1,1.02,0.94,1.21,0.95,1,1,1.06,1,1.04,0.94)
)

growth_df %>% group_by(shRNA,day) %>% summarise(mean = mean(val),
                                          sd = sd(val)) %>% data.frame() -> growth_df_plot
growth_df_plot$shRNA <- factor(as.character(growth_df_plot$shRNA), c("shluc", "sh1", "sh2"))

# make plot one
p2 <- ggplot(growth_df_plot, aes(x = day, y = mean, color = shRNA, group = shRNA)) +
  pretty_plot(fontsize = 8) +
  geom_line() +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  scale_color_manual(values = c("shluc" = "black", "sh2" = "orange3", "sh1" = "firebrick")) + L_border() +
  labs(x = "", y = "Fold Expansion") + theme(legend.position = "none")  +
  scale_y_continuous(expand = c(0,0), limits = c(0.8,2.5), breaks = c(1,1.5,2, 2.5)) +
  geom_hline(yintercept = 1, linetype = 2)

# 3) qPCR data
ann_df <- data.frame(
  shRNA = c(rep("shluc", 2), rep("sh1", 2),  rep("sh2", 2)),
  val = c(4.9,4.5,17.15,17.24,18.64,19)
)

ann_df %>% group_by(shRNA) %>% summarise(mean = mean(val),
                                          sd = sd(val)) %>% data.frame() -> ann_df_plot
ann_df_plot$shRNA <- factor(as.character(ann_df_plot$shRNA), c("shluc", "sh1", "sh2"))

# make plot one
p3 <- ggplot(ann_df_plot, aes(x = shRNA, y = mean, fill = shRNA)) +
  pretty_plot(fontsize = 8) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("shluc" = "black", "sh2" = "orange3", "sh1" = "firebrick")) + L_border() +
  labs(x = "", y = "AnnexinV+ %") + theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0))


cowplot::ggsave(cowplot::plot_grid(p1, p2, p3, nrow = 3),
                file = "threePanels.pdf", width = 1.3, height = 3.7)

