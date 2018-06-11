library(reshape2)
library(dplyr)
library(BuenColors)
library(stringr)
library(data.table)
  library(scales)

# Import RNA-seq
lapply(list.files("../../downloads/RNA_DESeq2"), function(x){
  tab <- data.frame(fread(paste0("../../downloads/RNA_DESeq2", "/", x), header = TRUE))
  s <- str_split_fixed(gsub("RNAseq-", "", x), "_", 3)
  pop1 <- s[1,1]
  pop2 <- s[1,2]
  count <- dim(tab)[1]
  data.frame(pop1 = pop1, pop2 = pop2, count = count)
}) %>% rbindlist() %>% data.frame() -> RNAseqCounts

# Import ATAC
lapply(list.files("../../downloads/ATAC_DESeq2"), function(x){
  tab <- data.frame(fread(paste0("../../downloads/ATAC_DESeq2", "/", x), header = TRUE))
  s <- str_split_fixed(gsub("ATACseq-", "", x), "_", 3)
  pop1 <- s[1,1]
  pop2 <- s[1,2]
  count <- dim(tab)[1]
  data.frame(pop1 = pop2, pop2 = pop1, count = -1*count)
}) %>% rbindlist() %>% data.frame() -> ATACseqCounts

# Combine and plot
total <- rbind(RNAseqCounts, ATACseqCounts)
total$pop1 <- factor(as.character(total$pop1),  levels = paste0("P", as.character(1:8)))
total$pop2 <- factor(as.character(total$pop2),  levels = rev(paste0("P", as.character(1:8))))

p1 <- ggplot(total, aes(x = pop1, y = pop2, fill = count)) + geom_tile() + pretty_plot() + 
  scale_fill_gradientn(colors = jdb_palette("solar_flare"), values =rescale(c(min(total$count), 0, max(total$count)))) +
  geom_abline(slope=-1, intercept = 9, linetype = 2)  + labs(x = "", y = "") + theme(legend.position = "none")

cowplot::ggsave(p1, file = "../../plots/DEseqHeatmap.pdf", width = 2.5, height = 2.5)

legendBottom <- g_legend( ggplot(total, aes(x = pop1, y = pop2, fill = count)) + geom_tile() + pretty_plot() + 
  scale_fill_gradientn(colors = rev(jdb_palette("solar_flare")[1:5]))+ theme(legend.position = "bottom"))
cowplot::ggsave(legendBottom, file = "bottom.pdf")
