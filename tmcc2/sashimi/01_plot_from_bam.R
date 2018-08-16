library(ggplot2)
library(reshape2)
library(Gviz)
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg19)
seqlevels(Hsapiens) <- gsub("chr", "", seqlevels(Hsapiens))
options(ucscChromosomeNames=FALSE)
options(R_MAX_NUM_DLLS=1000)

# Chromosome 6: 135,181,315-135,219,173 

PLOT_chr="1"
PLOT_min=205197038 - 50000
PLOT_max=205242471 + 50000
gtr <- GenomeAxisTrack()

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160",
                      "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")

P1 <- AlignmentsTrack("bams/aP1_TMCC2.bam", fill = "#3b82ae", col = "black",
                      name="P1", isPaired = TRUE, sashimiScore=20, chromosome=PLOT_chr, lwd.sashimiMax = 4)
P2 <- AlignmentsTrack("bams/aP2_TMCC2.bam", fill = "#547294", col = "black",
                      name="P2", isPaired = TRUE, sashimiScore=20, chromosome=PLOT_chr, lwd.sashimiMax = 4)
P3 <- AlignmentsTrack("bams/aP3_TMCC2.bam", fill = "#6d617a", col = "black",
                      name="P3", isPaired = TRUE, sashimiScore=20, chromosome=PLOT_chr, lwd.sashimiMax = 4)
P4 <- AlignmentsTrack("bams/aP4_TMCC2.bam", fill = "#865160", col = "black",
                      name="P4", isPaired = TRUE, sashimiScore=20, chromosome=PLOT_chr, lwd.sashimiMax = 4)
P5 <- AlignmentsTrack("bams/aP5_TMCC2.bam", fill = "#9f4046", col = "black",
                      name="P5", isPaired = TRUE, sashimiScore=20, chromosome=PLOT_chr, lwd.sashimiMax = 4)
P6 <- AlignmentsTrack("bams/aP6_TMCC2.bam", fill = "#b8302c", col = "black",
                      name="P6", isPaired = TRUE, sashimiScore=20, chromosome=PLOT_chr, lwd.sashimiMax = 4)
P7 <- AlignmentsTrack("bams/aP7_TMCC2.bam", fill = "#d11f12", col = "black",
                      name="P7", isPaired = TRUE, sashimiScore=20, chromosome=PLOT_chr, lwd.sashimiMax = 4)
P8 <- AlignmentsTrack("bams/aP8_TMCC2.bam", fill = "#de1705", col = "black",
                      name="P8", isPaired = TRUE, sashimiScore=20, chromosome=PLOT_chr, lwd.sashimiMax = 4)


if(TRUE){
  pdf(file = paste0("TMCC2_a.pdf"), width = 3, height = 3)
  plotTracks(trackList = list(P1, P2, P3, P4), start = c(PLOT_min),
             transcriptAnnotation= "symbol",from=PLOT_min, to=PLOT_max,
             type = c("coverage"), fontsize=6, col.id="black", col.axis="black", sizes= rep(1, 4),
             background.title = "white", innerMargin = 0, margin = 0)
  dev.off()
}

pdf(file = paste0("TMCC2_b.pdf"), width = 3, height = 3)
plotTracks(trackList = list(P5, P6), start = c(PLOT_min),
           transcriptAnnotation= "symbol",from=PLOT_min, to=PLOT_max,
           type = c("coverage"), fontsize=6, col.id="black", col.axis="black", sizes= rep(1, 2),
           background.title = "white", innerMargin = 0, margin = 0)
dev.off()

pdf(file = paste0("TMCC2_c.pdf"), width = 3, height = 3)
plotTracks(trackList = list(P7, P8), start = c(PLOT_min),
           transcriptAnnotation= "symbol",from=PLOT_min, to=PLOT_max,
           type = c("coverage"), fontsize=6, col.id="black", col.axis="black", sizes= rep(1, 2),
           background.title = "white", innerMargin = 0, margin = 0)
dev.off()
