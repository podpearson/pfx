# genotypeCallsFromGTas012.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


genotypeCallsFromGTas012 <- function(
  vcf,
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0)
) {
  GTs <- geno(vcf)[["GT"]]
  GTsInt <- matrix(GTsToIntMapping[GTs], nrow=nrow(GTs), dimnames=dimnames(GTs))
  return(GTsInt)
}
