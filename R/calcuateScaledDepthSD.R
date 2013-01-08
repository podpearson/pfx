# calcuateScaledDepthSD.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


calcuateScaledDepthSD <- function(
  vcf
) {
  RefReads <- matrix(
    sapply(geno(vcf)[["AD"]], function(x) x[1]),
    ncol=dim(geno(vcf)[["AD"]])[2],
    dimnames=dimnames(geno(vcf)[["AD"]])
  )
  FirstAltReads <- matrix(
    sapply(geno(vcf)[["AD"]], function(x) x[2]),
    ncol=dim(geno(vcf)[["AD"]])[2],
    dimnames=dimnames(geno(vcf)[["AD"]])
  )
  AllReads <- RefReads + FirstAltReads
  AllReads[is.na(AllReads)] <- 0
  scaledDepths <- scale(AllReads)
  scaledDepthSD <- apply(scaledDepths, 1, sd)
  return(scaledDepthSD)
}
