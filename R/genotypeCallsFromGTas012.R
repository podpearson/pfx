# genotypeCallsFromGTas012.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


genotypeCallsFromGTas012 <- function(
  vcf,
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0, "./."=0, "2"=0, "3"=0) # "./." is needed as sometimes this is output by GATK's UG (presumably a bug). "2", "3", needed for the case of multi-allelic sites
) {
  GTs <- geno(vcf)[["GT"]]
  GTsInt <- matrix(GTsToIntMapping[GTs], nrow=nrow(GTs), dimnames=dimnames(GTs))
  return(GTsInt)
}
