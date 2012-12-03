# convertGTsIntToParentBasedGTs.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


convertGTsIntToParentBasedGTs <- function(
  GTsToUse,
  IDparent1                   = dimnames(GTsInt)[[2]][1],
  IDparent2                   = dimnames(GTsInt)[[2]][2],
  reverseOrderOfSamples       = TRUE,
  return0asNA                 = FALSE
) {
  GTsFirstParent <- matrix(GTsToUse[, IDparent1], nrow=nrow(GTsToUse), ncol=ncol(GTsToUse))
  GTsSecondParent <- matrix(GTsToUse[, IDparent2], nrow=nrow(GTsToUse), ncol=ncol(GTsToUse))
  if(reverseOrderOfSamples) {
    GTsCFparents <- (
      ((GTsToUse == 0) * 0) +
      ((GTsToUse != 0 & GTsFirstParent !=0 & GTsSecondParent !=0 & GTsFirstParent != GTsSecondParent & GTsToUse == GTsFirstParent) * 1) +
      ((GTsToUse != 0 & GTsFirstParent !=0 & GTsSecondParent !=0 & GTsFirstParent != GTsSecondParent & GTsToUse == GTsSecondParent) * 2) +
      ((GTsToUse != 0 & (GTsFirstParent ==0 | GTsSecondParent ==0) & GTsFirstParent != GTsSecondParent & GTsToUse == GTsFirstParent) * 3) +
      ((GTsToUse != 0 & (GTsFirstParent ==0 | GTsSecondParent ==0) & GTsFirstParent != GTsSecondParent & GTsToUse == GTsSecondParent) * 4) +
      ((GTsToUse != 0 & GTsFirstParent !=0 & GTsSecondParent !=0 & GTsFirstParent == GTsSecondParent & GTsToUse == GTsFirstParent) * 5) +
      ((GTsToUse != 0 & GTsFirstParent !=0 & GTsSecondParent !=0 & GTsFirstParent == GTsSecondParent & GTsToUse != GTsFirstParent) * 6) +
      ((GTsToUse != 0 & GTsFirstParent ==0 & GTsSecondParent ==0 * 7))
    )[, seq(dim(GTsToUse)[2], 1, -1)] # reverse order so parents are last and hence appear at top of plots
  } else {
    GTsCFparents <- (
      ((GTsToUse == 0) * 0) +
      ((GTsToUse != 0 & GTsFirstParent !=0 & GTsSecondParent !=0 & GTsFirstParent != GTsSecondParent & GTsToUse == GTsFirstParent) * 1) +
      ((GTsToUse != 0 & GTsFirstParent !=0 & GTsSecondParent !=0 & GTsFirstParent != GTsSecondParent & GTsToUse == GTsSecondParent) * 2) +
      ((GTsToUse != 0 & GTsFirstParent !=0 & GTsSecondParent ==0 & GTsFirstParent != GTsSecondParent & GTsToUse == GTsFirstParent) * 3) +
      ((GTsToUse != 0 & GTsFirstParent ==0 & GTsSecondParent !=0 & GTsFirstParent != GTsSecondParent & GTsToUse == GTsSecondParent) * 4) +
      ((GTsToUse != 0 & GTsFirstParent !=0 & GTsSecondParent !=0 & GTsFirstParent == GTsSecondParent & GTsToUse == GTsFirstParent) * 5) +
      ((GTsToUse != 0 & GTsFirstParent !=0 & GTsSecondParent !=0 & GTsFirstParent == GTsSecondParent & GTsToUse != GTsFirstParent) * 6) +
      ((GTsToUse != 0 & GTsFirstParent ==0 & GTsSecondParent ==0 * 7))
    )
  }
  if(return0asNA) {
    GTsCFparents[GTsCFparents==0] <- NA
  }
  return(GTsCFparents)
}
