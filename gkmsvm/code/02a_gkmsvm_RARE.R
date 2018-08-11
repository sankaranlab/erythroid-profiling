library(seqinr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(gkmSVM)

if(FALSE){
  snps_df <- read.table("../../data/mendelianVariants.tsv", header = TRUE)
  pad <- 9
  
  # Make a GenomicRanges object of the regions of interest
  snps_df$start <- snps_df$BP - pad
  snps_df$end <- snps_df$BP + pad
  snps_gr <- makeGRangesFromDataFrame(snps_df, keep.extra.columns = TRUE)
  
  # Extract sequences from a reference genome
  sequence <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, snps_gr))
  
  # Annotate with the reference and the alternate alleles
  ref <-  paste0(substring(sequence, 1,9), snps_df$REF, substring(sequence, 11,19))
  alt <-  paste0(substring(sequence, 1,9), snps_df$ALT, substring(sequence, 11,19))
  names <- paste0(snps_df$CHR, ":", snps_df$BP, "_", snps_df$REF, "-", snps_df$ALT)
}

if(FALSE){
  write.fasta(as.list(ref), names, "../fasta/mendelianRare_ref.fasta",
              open = "w", nbchar = 60)
  
  write.fasta(as.list(alt), names, "../fasta/mendelianRare_alt.fasta",
              open = "w", nbchar = 60, as.string = FALSE)
}

# Function to compute delta 
do_gkmsvm_delta <- function(population){
  gkmsvm_delta('../fasta/mendelianRare_ref.fasta','../fasta/mendelianRare_alt.fasta',
               svmfnprfx=paste0('../kernel/', population),
               paste0('../variantPredictions/',population,'_mendelianRare.out'))
}

lapply(paste0("P", 1:8), do_gkmsvm_delta)


