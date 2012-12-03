# annotateSegregationStatus.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


annotateSegregationStatus <- function(
  vcf,
  parentalIDs                 = dimnames(vcf)[[2]][1:2]
) {
  GTsInt <- genotypeCallsFromGTas012(vcf)
  columnIndexesOfParents <- match(parentalIDs, dimnames(GTsInt)[[2]])
  columnIndexesOfProgeny <- setdiff(seq(along=dimnames(GTsInt)[[2]]), match(parentalIDs, dimnames(GTsInt)[[2]]))
  segregating <- (
    (GTsInt[, parentalIDs[1]] == 1 & GTsInt[, parentalIDs[2]] == 2) |
    (GTsInt[, parentalIDs[1]] == 2 & GTsInt[, parentalIDs[2]] == 1)
  )
  info(vcf) <- cbind(
    values(info(vcf)),
    DataFrame(
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
        DataFrame(Number="0", Type="Flag", Description="Is this a segregating site (i.e. do parents have different genotypes", row.names="SEGREGATING"),
        DataFrame(Number="1", Type="Integer", Description="Number of Mendelian errors (parents have same genotype, progeny has differenet genotype) in progeny", row.names="MendelianErrors")
      )
    )
  )
  return(vcf)
}
