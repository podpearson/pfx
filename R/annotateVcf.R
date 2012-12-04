# annotateVcf.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


annotateVcf <- function(
  vcf,
  parentalIDs                 = dimnames(vcf)[[2]][1:2]
) {
  GTsInt <- genotypeCallsFromGTas012(vcf)
  columnIndexesOfParents <- match(parentalIDs, dimnames(GTsInt)[[2]])
  columnIndexesOfProgeny <- setdiff(seq(along=dimnames(GTsInt)[[2]]), match(parentalIDs, dimnames(GTsInt)[[2]]))
#  segregating <- (
#    (GTsInt[, parentalIDs[1]] == 1 & GTsInt[, parentalIDs[2]] == 2) |
#    (GTsInt[, parentalIDs[1]] == 2 & GTsInt[, parentalIDs[2]] == 1)
#  )
  ADsArray <- array(
    unlist(geno(vcf)[["AD"]]),
    dim=c(2, dim(geno(vcf)[["AD"]])[1], dim(geno(vcf)[["AD"]])[2]),
    dimnames=list(c("Ref", "Nonref"), dimnames(geno(vcf)[["AD"]])[[1]], dimnames(geno(vcf)[["AD"]])[[2]])
#  Was originally doing this as below which is incorrect
#    dim=c(2, dim(geno(vcf)[["AD"]])[2], dim(geno(vcf)[["AD"]])[1]),
#    dimnames=list(c("Ref", "Nonref"), dimnames(geno(vcf)[["AD"]])[[2]], dimnames(geno(vcf)[["AD"]])[[1]])
  )
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
  MAF <- pmin(RefReads, FirstAltReads)/(RefReads+FirstAltReads)
#  MAF <- pmin(ADsArray[1,,], ADsArray[2,,])/(ADsArray[1,,]+ADsArray[2,,])
  meanVariantMAFs <- rowMeans(MAF, na.rm = TRUE)
  
  missingness <- apply(GTsInt, 1, function(x) length(which(x==0)))
  
  ADas0123 <- genotypeCallsFromADas0123(vcf)
  missingnessPerVariant <- apply(ADas0123, 1, function(x) length(which(is.na(x) | x==0)))
  heterozygosityPerVariant <- apply(ADas0123, 1, function(x) length(which(x==3)))
  
# debugging stuff - ignore
#  ADsArray[,1,2]
#  geno(vcf)[["AD"]][1,2]
#  ADsArray[,2,1]
#  geno(vcf)[["AD"]][2,1]
#  geno(vcf)[["AD"]][26,1]
#  geno(vcf)[["AD"]][27,1]
#  ADsArray[,26,1]
#  ADsArray[,27,1]
#  RefReads[26,1]
#  RefReads[27,1]
#  FirstAltReads[26,1]
#  FirstAltReads[27,1]
#  v2 <- unlist(geno(vcf)[["AD"]][2,])
#  ADsTemp <- array(
#    unlist(geno(vcf)[["AD"]][2, ]),
#    dim=c(2, 1, dim(geno(vcf)[["AD"]])[1]),
#    dimnames=list(c("Ref", "Nonref"), dimnames(geno(vcf)[["AD"]])[[2]][2], dimnames(geno(vcf)[["AD"]])[[1]])
#  )
#  browser()
  info(vcf) <- cbind(
    values(info(vcf)),
    DataFrame(
      meanMAF = meanVariantMAFs,
      missingness = missingness,
      missingness2 = missingnessPerVariant,
      heterozgosity = heterozygosityPerVariant,
      SEGREGATING=(
        (GTsInt[, parentalIDs[1]] == 1 & GTsInt[, parentalIDs[2]] == 2) |
        (GTsInt[, parentalIDs[1]] == 2 & GTsInt[, parentalIDs[2]] == 1)
      ),
      MendelianErrors=apply(
        GTsInt,
        1,
        function(genotypesForVariant) {
          sum(
            genotypesForVariant[parentalIDs[1]] == 1 &
              genotypesForVariant[parentalIDs[2]] == 1 &
              genotypesForVariant[columnIndexesOfProgeny] != 1
          ) +
          sum(
            genotypesForVariant[parentalIDs[1]] == 2 &
              genotypesForVariant[parentalIDs[2]] == 2 &
              genotypesForVariant[columnIndexesOfProgeny] != 2
          )
        }
      )
    )
  )
  exptData(vcf)[["header"]] <- VCFHeader(
    reference=reference(exptData(vcf)[["header"]]),
    samples=samples(exptData(vcf)[["header"]]),
    header=DataFrameList(
      META=meta(exptData(vcf)[["header"]]),
      FILTER=fixed(exptData(vcf)[["header"]])[["FILTER"]],
#      fixed=fixed(exptData(vcf)[["header"]]),
      FORMAT=geno(exptData(vcf)[["header"]]),
      INFO=rbind(
        info(exptData(vcf)[["header"]]),
        DataFrame(Number="1", Type="Float", Description="Mean across samples of proportion of minor allele reads (something like a heterozygosity score)", row.names="meanMAF"),
        DataFrame(Number="1", Type="Integer", Description="Number of samples with a missing genotype call", row.names="missingness"),
        DataFrame(Number="0", Type="Flag", Description="Is this a segregating site (i.e. do parents have different genotypes", row.names="SEGREGATING"),
        DataFrame(Number="1", Type="Integer", Description="Number of Mendelian errors (parents have same genotype, progeny has differenet genotype) in progeny", row.names="MendelianErrors")
      )
    )
  )
  return(vcf)
}
