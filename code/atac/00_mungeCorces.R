library(dplyr)
library(data.table)

# Now deal with expression
pp <- data.frame(fread("../../data/corces/atac_individualReplicates/panHeme_peaks_Preplicates.counts.tsv"))

# Import mapping / translate SRR files
mapdf <- read.table("../../data/corces/ATAC_Corces.txt", sep = "\t")
mapvec <- as.character(mapdf[,2]); names(mapvec) <- as.character(mapdf[,1])

# Import counts
pub <- data.frame(fread("../../data/corces/atac_individualReplicates/corces_replicates.counts.tsv"))
colnames(pub) <- gsub(".proatac", "", colnames(pub))

ATAC.counts.all <- cbind(pub, pp)
ATAC.celltype <- c(mapvec[colnames(pub)], stringr::str_split_fixed(colnames(pp), "_", 4)[1:28,3])

save(ATAC.counts.all, ATAC.celltype, file = "../../processed/serialized/ATACseq_wCorces.rda")
