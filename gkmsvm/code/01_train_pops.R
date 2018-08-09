library(gkmSVM)
suppressMessages(suppressWarnings(library(tools)))

do_gkmSVM <- function(pop){
  genNullSeqs('../data/topPeaks-P1.bed',nMaxTrials=10,xfold=1,
              genomeVersion='hg19',
              outputPosFastaFN=paste0('../fasta/',pop,'-positive.fa'),
              outputBedFN=paste0('../fasta/',pop,'-negative.bed'),
              outputNegFastaFN=paste0('../fasta/',pop,'-negative.fa'))
  
  gkmsvm_kernel(paste0('../fasta/',pop,'-positive.fa'),
                paste0('../fasta/',pop,'-negative.fa'),
                paste0('../kernel/', pop, ".kernel.out"))
  
  gkmsvm_trainCV(kernelfn = paste0('../kernel/', pop, ".kernel.out"),
                 posfn = paste0('../fasta/',pop,'-positive.fa'),
                 negfn = paste0('../fasta/',pop,'-negative.fa'), 
                 svmfnprfx=paste0('../kernel/', pop),
                 outputCVpredfn=paste0('../kernel/', pop, ".cvPred.out"),
                 outputROCfn=paste0('../kernel/', pop, ".roc.out"))
  
  gkmsvm_classify('../fasta/nr10mers.fa',
                  svmfnprfx=paste0('../kernel/', pop),
                  paste0('../kernel/', pop, ".weights.10mer.out"))
  pop
  
}

args <- commandArgs(trailingOnly = TRUE)
if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}

pops <- paste0("P", as.character(1:8))
idx <- as.numeric(args[i+1])

do_gkmSVM(pops[idx])
