library(GenomicAlignments)
library(rtracklayer)

makeBigwig <- function(pop){ 
  bamfile <- paste0("bams/a", pop, "_TMCC2.bam")
  alignment <- readGAlignments(bamfile,param = ScanBamParam(
    flag = scanBamFlag(), mapqFilter = 0))
  
  reads_coverage <- coverage(alignment)
  export.bw(reads_coverage, con = paste0("bigwigs/", pop, ".bw"))
}

lapply(paste0("P", as.character(1:8)), makeBigwig)