# convertGTsIntToParentBasedGTs.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


convertGTsIntToParentBasedGTs <- function(
  GTsToUse,
  IDparent1                   = dimnames(GTsToUse)[[2]][1],
  IDparent2                   = dimnames(GTsToUse)[[2]][2],
  reverseOrderOfSamples       = TRUE,
  return0asNA                 = FALSE,
  returnNonSegregatingAsNA    = FALSE
) {
  GTsFirstParent <- matrix(GTsToUse[, IDparent1], nrow=nrow(GTsToUse), ncol=ncol(GTsToUse))
  GTsSecondParent <- matrix(GTsToUse[, IDparent2], nrow=nrow(GTsToUse), ncol=ncol(GTsToUse))
  if(reverseOrderOfSamples) {
    GTsCFparents <- (
      ((GTsToUse == 0) * 0) +                                                                                                                   # mising genotype in sample
      ((GTsToUse != 0 & GTsFirstParent !=0 & GTsSecondParent !=0 & GTsFirstParent != GTsSecondParent & GTsToUse == GTsFirstParent) * 1) +       # segregating and same as first parent
      ((GTsToUse != 0 & GTsFirstParent !=0 & GTsSecondParent !=0 & GTsFirstParent != GTsSecondParent & GTsToUse == GTsSecondParent) * 2) +      # segregating and same as second parent
      ((GTsToUse != 0 & (GTsFirstParent ==0 | GTsSecondParent ==0) & GTsFirstParent != GTsSecondParent & GTsToUse == GTsFirstParent) * 3) +     # same as first parent but second parent missing
      ((GTsToUse != 0 & (GTsFirstParent ==0 | GTsSecondParent ==0) & GTsFirstParent != GTsSecondParent & GTsToUse == GTsSecondParent) * 4) +    # same as second parent but first parent missing
      ((GTsToUse != 0 & GTsFirstParent !=0 & GTsSecondParent !=0 & GTsFirstParent == GTsSecondParent & GTsToUse == GTsFirstParent) * 5) +       # non-segregating and same as parents
      ((GTsToUse != 0 & GTsFirstParent !=0 & GTsSecondParent !=0 & GTsFirstParent == GTsSecondParent & GTsToUse != GTsFirstParent) * 6) +       # non-segregating and different to parents (i.e. a Mendelian error)
      ((GTsToUse != 0 & GTsFirstParent ==0 & GTsSecondParent ==0 * 7))                                                                          # both parents missing
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
  if(returnNonSegregatingAsNA) {
    GTsCFparents[GTsCFparents>2] <- NA
  }
  return(GTsCFparents)
}
