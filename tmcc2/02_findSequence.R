library(GenomicRanges) 
library(BSgenome.Hsapiens.UCSC.hg19)
library(stringr)

# Make a data frame / GRanges of the gene region
gene_df <- data.frame(
 chr = c("chr1"),
 start = c(205185680),
 end = c(205253830),
 ID = c("TMCC2")
)

gene_gr <- makeGRangesFromDataFrame(gene_df, keep.extra.columns = TRUE)

# Extract sequences from the reference genome
sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gene_gr)

# String match -- returns coordinates
coords <- str_locate_all(as.character(sequence), fixed("AAGTGTTCGAGAAGA"))[[1]][1,] # unique string sequence from 364
location_shRNA_hg19 <- gene_df$start + coords
location_shRNA_hg19