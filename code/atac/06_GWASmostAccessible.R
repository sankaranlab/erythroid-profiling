library(dplyr)
library(GenomicRanges)
library(BiocParallel)
library(data.table)
library(ComplexHeatmap)
library(BuenColors)
set.seed(14651)
register(MulticoreParam(2))

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")
color_map_all <- c(ejc_color_maps, eryth_color_maps)

# Import
bed <- data.frame(fread("../../data/ATAC_data/ery_only.bed", header = FALSE))
colnames(bed) <- c("chr", "start", "end")
peak_gr <- makeGRangesFromDataFrame(bed)

# Import raw expression values
pops <- c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8)))
total <-  data.matrix(data.frame(fread("../../data/corces/atac_combinedPopulations/panHeme.counts.tsv")))[,pops]

# Normalize
cpm <- log2(sweep(total, 2, colSums(total), FUN="/") * 1000000 + 1)

# Find overlaps with SNPS
import_to_GR <- function(trait){
  file <- paste0("../../data/FMsnps/", trait, "_PP001_betas.bed")
  tab <- read.table(file)
  tab <- tab[tab$V5 > 0.50,]
  colnames(tab) <- c("chr", "start", "end", "region", "PP", "Beta", "SE", "Z")
  gr <- makeGRangesFromDataFrame(tab[,c(1,2,3,5,8)], keep.extra.columns = TRUE)
  ov1 <- findOverlaps(gr, peak_gr)
  terminal_df <- data.frame(gr[queryHits(ov1)], round(cpm[subjectHits(ov1),], 3))
  terminal_df$trait <- trait
  terminal_df
}

lapply(c("HCT", "HGB", "MCH", "MCHC", "MEAN_RETIC_VOL", "RBC_COUNT", "RETIC_COUNT"), import_to_GR) %>% 
  rbindlist() %>% data.frame() -> allDF

countsDFaccessiblility <- unique(allDF[,c("seqnames", "start", "end", pops)])
peaksKeep <- paste0(countsDFaccessiblility$seqnames, ":",
                    as.character(countsDFaccessiblility$start), "-",
                    as.character(countsDFaccessiblility$end))


# Z score
log2cpm <- countsDFaccessiblility[,4:15]
log2cpm_z <- sapply(1:dim(log2cpm)[1],function(i){
  x <- data.matrix(log2cpm)[i,]
  (x - mean(x))/sd(x)
}) %>% t()

# Max - min
normalized <- lapply(1:dim(log2cpm)[1], function(i){
  (log2cpm[i,] - min(log2cpm[i,]))/(max(log2cpm[i,])- min(log2cpm[i,]))
}) %>% rbindlist() %>% data.matrix() 

km <- kmeans(normalized, centers = 3, nstart = 10000)
km.cluster <- factor(km$cluster)

hm <- Heatmap(normalized, col=as.character(jdb_palette("brewer_spectra",type="continuous")),
        cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 0),
        column_names_gp = gpar(fontsize = 6),
        split = km.cluster, show_heatmap_legend = FALSE,
        name = "Peak\nAccessibility")
#hm

# Find most accessible
sapply(1:dim(log2cpm)[1], function(i){
  colnames(log2cpm)[which.max(log2cpm[i,])]
}) -> popVec

lapply(pops, function(Cell){
  data.frame(Cell = Cell, Count = sum(popVec == Cell))
}) %>% rbindlist() %>% data.frame() -> mdf

p1 <- ggplot(mdf, aes(x = Cell, y = Count, fill = Cell)) +
  geom_bar(stat = "identity", color = "black") + pretty_plot(fontsize = 8) + L_border() + 
   labs(x = "", y = "# SNPs Most Accessible", color = "") + theme(legend.position = "none") +
  scale_fill_manual(values = color_map_all)

ggsave(p1, filename = "../../plots/mostAccessible.pdf", width = 3.5, height = 1.5)

