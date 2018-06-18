library(dplyr)
library(irlba)
library(annotables)
library(BuenColors)
library(chromVAR)
library(chromVARmotifs)
library(SummarizedExperiment)
library(BiocParallel)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(circlize)
library(ComplexHeatmap)
register(MulticoreParam(2))

set.seed(14651)

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")
all_color_maps <- c(ejc_color_maps[c("HSC", "MPP", "CMP")], "MEP" = "#FF81AF" ,eryth_color_maps)

# Import raw expression values
load("../../processed/serialized/ATACseq_wCorces.rda")

# chromVAR Setup
ATAC.counts <- data.matrix(data.frame(data.table::fread("../../data/corces/atac_combinedPopulations/panHeme.counts.tsv")))
pops.ordered <- c("HSC", "MPP", "CMP", "MEP", paste0("P", as.character(1:8)))

SE <- SummarizedExperiment(
 rowRanges =  diffloop::bedToGRanges("../../data/corces/panHeme.bed"),
 colData = data.frame(Sample =pops.ordered, population = pops.ordered),
 assays = list(counts = Matrix::Matrix(ATAC.counts[,pops.ordered]))
)
SE <- filterPeaks(SE)

SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
data("human_pwms_v2")

# Match and deviate
motif_ix <- matchMotifs(human_pwms_v2, SE, genome = BSgenome.Hsapiens.UCSC.hg19)
names(assays(motif_ix)) <- "weights"
dev <- gchromVAR::computeWeightedDeviations(object = SE, weights = motif_ix)

# Do Variability
bagged <- chromVARxx::bagDeviations(dev, 0.8, "human")
variabilityAll <- computeVariability(bagged)

data.frame(variabilityAll) %>% arrange(desc(variability))-> motifVariability
motifVariability$rank <- 1:dim(variabilityAll)[1]
p1 <- ggplot(motifVariability, aes(x = rank, y = variability)) + geom_point(size = 0.75) +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "Motif Rank", y = "Accessibility Variability")
cowplot::ggsave(p1, file = "../../plots/motifVariability.pdf", width = 2.5, height = 2.5)

m <- bagged[rownames(variabilityAll)[variabilityAll$variability>9.9][1:20],]@assays$data$z
rownames(m) <- as.character(variabilityAll$name[variabilityAll$variability>9.9][1:20])
m <- ifelse(m > 20, 20, m)
m <- ifelse(m < -20, -20, m)


ha_col <- HeatmapAnnotation(cell = as.character(colData(bagged)$population),
                            col = list(cell = all_color_maps))
pdf("../../plots/deviations_wCorces.pdf", width = 6, height = 2.5)
par(cex.main=0.8,mar=c(2,2,2,2))

hm <- Heatmap(m, col=as.character(jdb_palette("solar_extra",type="continuous")),
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              show_row_names = TRUE, 
              cluster_columns = FALSE,
              show_column_names = FALSE)
hm
dev.off()
