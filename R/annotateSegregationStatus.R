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
  columnIndexesOfProgeny <- setdiff(seq(along=dimnames(GTsInt)[[2]], match(parentalIDs, dimnames(GTsInt)[[2]]))
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
              genotypesForVariant[genotypesForVariant] == 1 &
              genotypesForVariant[columnIndexesOfProgeny] != 1
          ) +
          sum(
            genotypesForVariant[parentalIDs[1]] == 2 &
              genotypesForVariant[genotypesForVariant] == 2 &
              genotypesForVariant[columnIndexesOfProgeny] != 2
          )
        }
      )
    )
  )
  browser()
}
