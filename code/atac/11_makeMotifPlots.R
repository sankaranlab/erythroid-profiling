library(ggseqlogo)
library(ggplot2)
library(chromVAR)
library(chromVARmotifs)

data("human_pwms_v1")
motifs <- chromVAR::getJasparMotifs()

#BCL11A
bcl11a <- human_pwms_v1[["ENSG00000119866_LINE848_BCL11A_D"]]
p_motif <- ggplot() + geom_logo( exp(bcl11a@profileMatrix) * 0.25 ) + theme_logo()
cowplot::ggsave(p_motif, file = "../../plots/motifs/BCL11a.pdf", width = 5, height = 2)

#TAL1
tal1 <- human_pwms_v1[["ENSG00000162367_LINE243_TAL1_D_N1"]]
p_motif <- ggplot() + geom_logo( exp(tal1@profileMatrix) * 0.25 ) + theme_logo()
cowplot::ggsave(p_motif, file = "../../plots/motifs/TAL1.pdf", width = 5, height = 2)

#GATA1
gata1 <-  human_pwms_v1[["ENSG00000102145_LINE2073_GATA1_D_N1"]]
p_motif <- ggplot() + geom_logo( exp(gata1@profileMatrix) * 0.25 ) + theme_logo()
cowplot::ggsave(p_motif, file = "../../plots/motifs/GATA1.pdf", width = 5, height = 2)


#CTCF from JASPAR
ctcf <- motifs$MA0139.1_CTCF
p_motif <- ggplot() + geom_logo( ctcf@profileMatrix ) + theme_logo()
cowplot::ggsave(p_motif, file = "../../plots/motifs/CTCF.pdf", width = 5, height = 1.5)


#GATA1-TAL1 from JASPAR
gt <- motifs$`MA0140.2_GATA1::TAL1`
p_motif <- ggplot() + geom_logo( gt@profileMatrix ) + theme_logo()
cowplot::ggsave(p_motif, file = "../../plots/motifs/GATA1-TAL1.pdf", width = 5, height = 1.5)


