library(GenomicRanges) 
library(BSgenome.Hsapiens.UCSC.hg19)
library(stringr)

# Make a data frame / GRanges of the gene region
gene_df <- data.frame(
 chr = c("chr1"),
 start = c(205225000),
 end =   c(205225300),
 ID = c("TMCC2")
)

gene_gr <- makeGRangesFromDataFrame(gene_df, keep.extra.columns = TRUE)

# Extract sequences from the reference genome
sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gene_gr)
as.character(sequence)