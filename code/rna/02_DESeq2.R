library(dplyr)
library(DESeq2)
library(BiocParallel)
library(IHW)
library(annotables)
register(MulticoreParam(2))

# Filter sex chromosomes
mdf <- merge(read.table("../../data/genes.tsv", header = FALSE), grch37, by.x = "V1", by.y = "ensgene")
symbol_keep <- unique(mdf[mdf$chr %in% c(as.character(1:22), "X"),"V2"])

# Import raw expression values
raw <- read.table("../../processed/RNAseq_rawGeneCounts.tsv", header = TRUE)
raw <- raw[raw[,29] %in% symbol_keep,]
RNA.counts <- raw[,1:28]
meta <- stringr::str_split_fixed(colnames(RNA.counts), "_", 4)


# RNA DEseq2 setup
RNA.counts.df <- as.data.frame(RNA.counts)

# Establish column data
RNA.condition <- meta[,3]
colData <- as.data.frame(RNA.condition)
row.names(colData) <- colnames(RNA.counts.df)

# Run DEseq2
RNA.dds <- DESeqDataSetFromMatrix(countData = RNA.counts.df, colData = colData, design = ~ RNA.condition)
RNA.dds <- DESeq(RNA.dds, parallel = TRUE)

# All pairwise comparisons
comp <- expand.grid(paste0("P", seq(1:8)), paste0("P", seq(1:8)), stringsAsFactors = FALSE)
comp <- subset(comp, Var1!=Var2)
comp$cond <- "RNA.condition"
comp <- comp %>% 
  mutate(key = paste0(pmin(Var1, Var2), "_" ,pmax(Var1, Var2), sep = "")) %>%
  dplyr::distinct(key, .keep_all=TRUE)

# Function to do all pair-wise DEseq comparisons
sig.pairwise <- function(x) {
  print(x)
  # Do the DESeq contrast
  resSig.RNA <- data.frame(results(RNA.dds, contrast=c(comp[x,3], comp[x,2], comp[x,1]), parallel = TRUE, filterFun=ihw, alpha = 0.01))
  resSig.RNA$gene <- raw[,29]
  resSig.RNA %>% filter(complete.cases(.)) %>% arrange(padj) %>% filter(padj < 0.01) -> out
  
  # Do some rounding
  out$baseMean <- round( out$baseMean, 1)
  out$log2FoldChange <- round( out$log2FoldChange, 1)
  out$pvalue <-sprintf("%.3e", out$pvalue)
  out$padj <-sprintf("%.3e", out$padj)
   
  # Output table
  write.table(out[,c("gene", "baseMean", "log2FoldChange", "pvalue", "padj")],
              file = paste0("../../downloads/RNA_DESeq2/RNAseq-", comp[x,4], "_Padj01.tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  comp[x,4]
  
}

all.pairwise.P1 <- function(x) {
  print(x)
  # Do the DESeq contrast
  resSig.RNA <- data.frame(results(RNA.dds, contrast=c(comp[x,3], comp[x,2], comp[x,1]), parallel = TRUE, filterFun=ihw, alpha = 0.01))
  resSig.RNA$gene <- raw[,29]
  resSig.RNA %>% filter(complete.cases(.)) %>% arrange(padj) -> out
  
  # Do some rounding
  out$baseMean <- round( out$baseMean, 1)
  out$log2FoldChange <- round( out$log2FoldChange, 1)
  out$pvalue <-sprintf("%.3e", out$pvalue)
  out$padj <-sprintf("%.3e", out$padj)
   
  # Output table
  write.table(out[,c("gene", "baseMean", "log2FoldChange", "pvalue", "padj")],
              file = paste0("../../processed/all_DEseq2/", comp[x,4], "_deseq2full.tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  comp[x,4]
  
}

dumb <- lapply(1:7, all.pairwise.P1)
dumb <- lapply(which(comp$Var1 == "P8"), all.pairwise.P1)


dim(comp)
dumb <- lapply(1:dim(comp)[1], sig.pairwise)
