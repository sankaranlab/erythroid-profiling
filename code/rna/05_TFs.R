library(data.table)
library(dplyr)
library(BuenColors)
library(plotly)

lapply(list.files("../../processed/all_DEseq2/"), function(x){
  df <- read.table(paste0("../../processed/all_DEseq2/", x), header = TRUE)
  df$timepoint <- stringr::str_split_fixed(x, "_", 3)[1,2]
  df <- df[df$baseMean > 500,]
  df[,c("gene", "log2FoldChange", "timepoint")]
}) %>% rbindlist() %>% data.frame() -> allFC

# Import human transcription factors
tf_df <-data.frame( fread("../../data/humanTFs.txt", sep = "\t", header = TRUE))
transcription_factors <- as.character(tf_df[,2])

go <- allFC[as.character(allFC$gene) %in% transcription_factors, ]
go <- rbind(go, data.frame(gene = unique(as.character(go$gene)), log2FoldChange = 0, timepoint = "P1"))
go$gene <- factor(as.character(go$gene), levels = unique(as.character(go$gene)))

p1 <- ggplot(go, aes(x = timepoint, y = log2FoldChange*-1, group = gene)) +
  pretty_plot(fontsize = 8) +
  geom_line(size = 0.25)  + L_border() +
  labs(x = "", y = "log2FC over P1") + theme(legend.position = "none")

plotly_plot <- ggplotly(p1) 
htmlwidgets::saveWidget(plotly::as_widget(plotly_plot), "TF-geneExpression.html")


