# genotypeCallsFromADas0123.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


genotypeCallsFromADas0123 <- function(
  vcf
) {
  ADsArray <- array(
    unlist(geno(vcf)[["AD"]]),
    dim=c(2, dim(geno(vcf)[["AD"]])[2], dim(geno(vcf)[["AD"]])[1]),
    dimnames=list(c("Ref", "Nonref"), dimnames(geno(vcf)[["AD"]])[[2]], dimnames(geno(vcf)[["AD"]])[[1]])
  )
  ADas0123 <- t(
    (((ADsArray[1, , ] + ADsArray[2, , ]) < 5) * 0) +
    (((ADsArray[1, , ] + ADsArray[2, , ]) >= 5 & ADsArray[2, , ] <=1) * 1) +
    (((ADsArray[1, , ] + ADsArray[2, , ]) >= 5 & ADsArray[1, , ] <=1) * 2) +
    (((ADsArray[1, , ] + ADsArray[2, , ]) >= 5 & ADsArray[1, , ] >=2 & ADsArray[2, , ] >=2) * 3)
  )
  return(ADas0123)
}
