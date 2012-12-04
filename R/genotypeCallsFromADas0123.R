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
  ADas0123 <-
    (((RefReads + FirstAltReads) < 5) * 0) +
    (((RefReads + FirstAltReads) >= 5 & FirstAltReads <=1) * 1) +
    (((RefReads + FirstAltReads) >= 5 & RefReads <=1) * 2) +
    (((RefReads + FirstAltReads) >= 5 & RefReads >=2 & FirstAltReads >=2) * 3)
  

#  ADsArray <- array(
#    unlist(geno(vcf)[["AD"]]),
#    dim=c(2, dim(geno(vcf)[["AD"]])[1], dim(geno(vcf)[["AD"]])[2]),
#    dimnames=list(c("Ref", "Nonref"), dimnames(geno(vcf)[["AD"]])[[1]], dimnames(geno(vcf)[["AD"]])[[2]])
#  )
#  ADsArray <- array(
#    unlist(geno(vcf)[["AD"]]),
#    dim=c(2, dim(geno(vcf)[["AD"]])[2], dim(geno(vcf)[["AD"]])[1]),
#    dimnames=list(c("Ref", "Nonref"), dimnames(geno(vcf)[["AD"]])[[2]], dimnames(geno(vcf)[["AD"]])[[1]])
#  )
#  ADas0123 <-
#    (((ADsArray[1, , ] + ADsArray[2, , ]) < 5) * 0) +
#    (((ADsArray[1, , ] + ADsArray[2, , ]) >= 5 & ADsArray[2, , ] <=1) * 1) +
#    (((ADsArray[1, , ] + ADsArray[2, , ]) >= 5 & ADsArray[1, , ] <=1) * 2) +
#    (((ADsArray[1, , ] + ADsArray[2, , ]) >= 5 & ADsArray[1, , ] >=2 & ADsArray[2, , ] >=2) * 3)
#  ADas0123 <- t(
#    (((ADsArray[1, , ] + ADsArray[2, , ]) < 5) * 0) +
#    (((ADsArray[1, , ] + ADsArray[2, , ]) >= 5 & ADsArray[2, , ] <=1) * 1) +
#    (((ADsArray[1, , ] + ADsArray[2, , ]) >= 5 & ADsArray[1, , ] <=1) * 2) +
#    (((ADsArray[1, , ] + ADsArray[2, , ]) >= 5 & ADsArray[1, , ] >=2 & ADsArray[2, , ] >=2) * 3)
#  )
  return(ADas0123)
}
