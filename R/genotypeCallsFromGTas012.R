# genotypeCallsFromGTas012.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


genotypeCallsFromGTas012 <- function(
  vcf,
  GTsToIntMapping             = c("0"=1, "0/0"=1, "0|0"=1, "1"=2, "1/1"=2, "1|1"=2, "."=0, "./."=0, "./."=0, "2"=0, "3"=0, "0/1"=0, "1/0"=0, "0|1"=0, "1|0"=0), # "./." is needed as sometimes this is output by GATK's UG (presumably a bug). "2", "3", needed for the case of multi-allelic sites
  replaceRownamesWithChromPos = TRUE
) {
  GTs <- geno(vcf)[["GT"]]
  if(replaceRownamesWithChromPos) {
    GTsInt <- matrix(
      GTsToIntMapping[GTs],
      nrow=nrow(GTs),
      dimnames=list(
        paste(as.character(seqnames(vcf)), start(rowData(vcf)), sep=":"),
        dimnames(GTs)[[2]]
      )
    )
  } else {
    GTsInt <- matrix(
      GTsToIntMapping[GTs],
      nrow=nrow(GTs),
      dimnames=dimnames(GTs)
    )
  }
  return(GTsInt)
}
