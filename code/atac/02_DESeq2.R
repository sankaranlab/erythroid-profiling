library(dplyr)
library(DESeq2)
library(BiocParallel)
library(IHW)
library(annotables)
register(MulticoreParam(2))

# Filter sex chromosomes
bed <- read.table("../../data/ATAC_data/ery_only.bed", header = FALSE)
  
  
# Import raw expression values
raw <- read.table("../../data/ATAC_data/ery_only.counts.tsv", header = TRUE)
ATAC.counts <- data.matrix(raw[,1:28])
meta <- stringr::str_split_fixed(colnames(ATAC.counts), "_", 4)

# RNA DEseq2 setup
counts.df <- as.data.frame(ATAC.counts)

# Establish column data
ATAC.condition <- meta[,3]
colData <- as.data.frame(ATAC.condition)
row.names(colData) <- colnames(counts.df)

# Run DEseq2
RNA.dds <- DESeqDataSetFromMatrix(countData = counts.df, colData = colData, design = ~ ATAC.condition)
RNA.dds <- DESeq(RNA.dds, parallel = TRUE)

# All pairwise comparisons
comp <- expand.grid(paste0("P", seq(1:8)), paste0("P", seq(1:8)), stringsAsFactors = FALSE)
comp <- subset(comp, Var1!=Var2)
comp$cond <- "ATAC.condition"
comp <- comp %>% 
  mutate(key = paste0(pmin(Var1, Var2), "_" ,pmax(Var1, Var2), sep = "")) %>%
  dplyr::distinct(key, .keep_all=TRUE)

# Function to do all pair-wise DEseq comparisons
sig.pairwise <- function(x) {
  print(x)
  # Do the DESeq contrast
  resSig.ATAC <- data.frame(results(RNA.dds, contrast=c(comp[x,3], comp[x,2], comp[x,1]), parallel = TRUE, filterFun=ihw, alpha = 0.01))
  resSig.ATAC$chr <- bed[,1]
  resSig.ATAC$start <- bed[,2]
  resSig.ATAC$end <- bed[,3]
  
  resSig.ATAC %>% filter(complete.cases(.)) %>% arrange(padj) %>% filter(padj < 0.01) -> out
  
  # Do some rounding
  out$baseMean <- round( out$baseMean, 1)
  out$log2FoldChange <- round( out$log2FoldChange, 1)
  out$pvalue <-sprintf("%.3e", out$pvalue)
  out$padj <-sprintf("%.3e", out$padj)
   
  # Output table
  write.table(out[,c("chr", "start", "end", "baseMean", "log2FoldChange", "pvalue", "padj")],
              file = paste0("../../downloads/ATAC_DESeq2/ATACseq-", comp[x,4], "_Padj01.tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  comp[x,4]
  
}

dim(comp)
dumb <- lapply(1:dim(comp)[1], sig.pairwise)
