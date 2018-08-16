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
PLOT_min=205195038
PLOT_max=205244471
gtr <- GenomeAxisTrack()

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160",
                      "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")

P1 <- DataTrack("bigwigs/P1.bw", fill = "#3b82ae", col = "#3b82ae",  type = "h",
                name="P1",  chromosome=PLOT_chr, ylim = c(0,10000))
P2 <- DataTrack("bigwigs/P2.bw", fill = "#547294", col = "#547294", type = "h",
                name="P2",  chromosome=PLOT_chr, ylim = c(0,10000))
P3 <- DataTrack("bigwigs/P3.bw", fill = "#6d617a", col = "#6d617a", type = "h",
                name="P3",  chromosome=PLOT_chr, ylim = c(0,10000))
P4 <- DataTrack("bigwigs/P4.bw", fill = "#865160", col = "#865160", type = "h",
                name="P4",  chromosome=PLOT_chr, ylim = c(0,10000))
P5 <- DataTrack("bigwigs/P5.bw", fill = "#9f4046", col = "#9f4046", type = "h",
                name="P5",  chromosome=PLOT_chr, ylim = c(0,10000))
P6 <- DataTrack("bigwigs/P6.bw", fill = "#b8302c", col = "#b8302c", type = "h",
                name="P6",  chromosome=PLOT_chr, ylim = c(0,10000))
P7 <- DataTrack("bigwigs/P7.bw", fill = "#d11f12", col = "#d11f12", type = "h",
                name="P7",  chromosome=PLOT_chr, ylim = c(0,10000))
P8 <- DataTrack("bigwigs/P8.bw", fill = "#de1705", col = "#de1705", type = "h",
                name="P8",  chromosome=PLOT_chr, ylim = c(0,10000))

pdf(file = paste0("TMCC2_a.pdf"), width = 3.5, height = 1.7)
plotTracks(trackList = list(P1, P2,P3, P4, P5, P6, P7, P8), start = c(PLOT_min), #, 
           from=PLOT_min, to=PLOT_max,
           fontsize=6, col.id="black", col.axis="black", sizes= rep(1, 8),
           background.title = "white", innerMargin = 0, margin = 0)
dev.off()
