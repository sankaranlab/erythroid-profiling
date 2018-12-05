library(dplyr)
library(DESeq2)
library(BiocParallel)
library(IHW)
library(annotables)
register(MulticoreParam(2))

# Import raw expression values
raw <- read.table("output/TMCC2-RNAseq_rawGeneCounts.tsv", header = TRUE)
RNA.counts <- raw[,1:7]
genes <- raw[,8]
meta <- data.frame(gRNA = c("NC", "NC", "NC", rep("g-pooled", 4)))

# RNA DEseq2 setup
RNA.counts.df <- as.data.frame(RNA.counts)

# Establish column data
RNA.condition <- meta[,1]
colData <- as.data.frame(RNA.condition)
row.names(colData) <- colnames(RNA.counts.df)

# Run DEseq2
RNA.dds <- DESeqDataSetFromMatrix(countData = RNA.counts.df, colData = colData, design = ~ RNA.condition)
RNA.dds <- DESeq(RNA.dds, parallel = TRUE)

# All pairwise comparisons
comp <- data.frame(
  Var1 = c("g-pooled"),
  Var2 = c("NC"), stringsAsFactors = FALSE
)
comp$cond <- "RNA.condition"
comp <- comp %>% 
  mutate(key = paste0(Var1, "_", Var2, sep = "")) %>%
  dplyr::distinct(key, .keep_all=TRUE)

# Function to do all pair-wise DEseq comparisons
sig.pairwise <- function(x) {
  print(x)
  # Do the DESeq contrast
  resSig.RNA <- data.frame(results(RNA.dds, contrast=c(comp[x,3], comp[x,2], comp[x,1]), parallel = TRUE, filterFun=ihw, alpha = 0.01))
  resSig.RNA$gene <- genes
  resSig.RNA %>% filter(complete.cases(.)) %>% arrange(padj) -> out
  
  # Do some rounding
  out$baseMean <- round( out$baseMean, 1)
  out$log2FoldChange <- round( out$log2FoldChange, 4)
  out$pvalue <-sprintf("%.9e", out$pvalue)
  out$padj <-sprintf("%.3e", out$padj)
  
  # Output table
  write.table(out[,c("gene", "baseMean", "log2FoldChange", "pvalue", "padj")],
              file = paste0("output/RNAseq-", comp[x,4], "_Padj01.tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  comp[x,4]
  
}

dim(comp)
dumb <- lapply(1:dim(comp)[1], sig.pairwise)
