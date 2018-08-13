library(ggplot2)
library(reshape2)
library(Gviz)
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg19)
#seqlevels(Hsapiens) <- gsub("chr", "", seqlevels(Hsapiens))
options(ucscChromosomeNames=TRUE)
options(R_MAX_NUM_DLLS=1000)

# Chromosome 6: 135,181,315-135,219,173 

PLOT_chr="chr6"
PLOT_min=135515056-500
PLOT_max=135515494 +500
gtr <- GenomeAxisTrack()
biomTrack <- Gviz::BiomartGeneRegionTrack(genome = "hg19", chromosome = PLOT_chr,
                                          start = PLOT_min, end = PLOT_max, name = "UCSC")
dBRD9 <- AlignmentsTrack("dBRD9_6d.bam",
                              name="dBRD9_6D", isPaired = FALSE, sashimiScore=20, chromosome=PLOT_chr, lwd.sashimiMax = 4)
DMSO <- AlignmentsTrack("DMSO_6d.bam",
                            name="DMSO_6D", isPaired = FALSE, sashimiScore=20, chromosome=PLOT_chr, lwd.sashimiMax = 4)

pdf(file = paste0("MYB.pdf"), width = 10, height = 10)
plotTracks(trackList = list(gtr, biomTrack, dBRD9, DMSO), start = c(PLOT_min),
           collapseTranscripts="all", transcriptAnnotation= "symbol",from=PLOT_min, to=PLOT_max,
           type = c("coverage", "sashimi"), fontsize=10, col.id="black", col.axis="black", sizes=c(0.1,0.2, 0.35,0.35),
           background.title = "white", innerMargin = 0, margin = 0)
dev.off()
