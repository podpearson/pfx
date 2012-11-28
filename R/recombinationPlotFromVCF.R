# recombinationPlotFromVCF.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


recombinationPlotFromVCF <- function(
  vcfFilename                 = "/data/galton/users/rpearson/crossesTesting/release/7g8xGb4-qcPlusSamples-0.1.vcf.gz",
  chromosome                  = "MAL4",
  parentalIDs                 = NULL, # this default will use IDs of first two samples in VCF. This can be over-ridden by supplying IDs, e.g. recombinationPlotFromVCF(parentalIDs=c("ERR027099", "ERR027100"))
  shouldRemoveInvariant       = TRUE,
  removeMissingInParents      = c("either", "both", "ignore"),
  shouldRemoveMendelianErrors = FALSE,
  filtersToRemove             = NULL,
  samplesToRemove             = NULL,
  ...
) {
  GTsInt <- genotypeCallsFromGTas012(
    createSingleChromosomeVariantSitesRdaFile(
      vcfFilename                 = vcfFilename,
      chromosome                  = chromosome,
      shouldRemoveInvariant       = shouldRemoveInvariant,
      filtersToRemove             = filtersToRemove,
      samplesToRemove             = samplesToRemove,
      overwriteExisting           = TRUE,
      saveAsRobjectFile           = FALSE
    )
  )
  if(is.null(parentalIDs)) {
    parentalIDs <- dimnames(GTsInt)[[2]][1:2]
  }
  if(!identical(parentalIDs, dimnames(GTsInt)[[2]][1:2])) {
    GTsInt <- cbind(
      GTsInt[, parentalIDs],
      GTsInt[, -which(dimnames(GTsInt)[[2]] %in% parentalIDs)]
    )
  }
  if(removeMissingInParents[1] == "either") {
    if(shouldRemoveMendelianErrors) {
      rowsToUse <- (!GTsInt[, 1]==0 & !GTsInt[, 2]==0) & !(GTsInt[, 1] == GTsInt[, 2])
    } else {
      rowsToUse <- !GTsInt[, 1]==0 & !GTsInt[, 2]==0
    }
    GTsInt <- GTsInt[rowsToUse, ]
  }
  if(removeMissingInParents[1] == "both") {
    if(shouldRemoveMendelianErrors) {
      rowsToUse <- (!GTsInt[, 1]==0 | !GTsInt[, 2]==0) & !(GTsInt[, 1] == GTsInt[, 2])
    } else {
      rowsToUse <- !GTsInt[, 1]==0 | !GTsInt[, 2]==0
    }
    GTsInt <- GTsInt[rowsToUse, ]
  }

  recombinationPlot(
    convertGTsIntToParentBasedGTs(
      GTsInt
    ),
    ...
  )
}
