# recombinationPlotFromVCF.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################

# pdf("MAL4_7g8xGb4_release0.1.pdf", width=14, height=4)
# recombinationPlotFromVCF()
# dev.off()

recombinationPlotFromVCF <- function(
  vcfFilename                 = "/data/galton/users/rpearson/crossesTesting/release/7g8xGb4-qcPlusSamples-0.1.vcf.gz",
  chromosome                  = "MAL4",
  parentalIDs                 = NULL, # this default will use IDs of first two samples in VCF. This can be over-ridden by supplying IDs, e.g. recombinationPlotFromVCF(parentalIDs=c("ERR027099", "ERR027100"))
  shouldRemoveInvariant       = TRUE,
  removeMissingInParents      = c("either", "both", "ignore"),
  shouldRemoveMendelianErrors = TRUE,
  regionsMask                 = varRegions_v2(), # will remove any variants in these regions. Set to NULL if you don't want to mask any variants out in this way
  filtersToRemove             = c("Heterozygous", "Missingness"),     # set to NULL if you want to keep all variants
  samplesToRemove             = NULL,
  verbose                     = TRUE,
  ...
) {
  require(VariantAnnotation)
  if(verbose) {
    cat("recombinationPlotFromVCF: reading data from ", vcfFilename, "\n", sep="")
  }
  GTsInt <- genotypeCallsFromGTas012(
    createSingleChromosomeVariantSitesRdaFile(
      vcfFilename                 = vcfFilename,
      chromosome                  = chromosome,
      shouldRemoveInvariant       = shouldRemoveInvariant,
      regionsMask                 = regionsMask,
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
    if(verbose) {
      cat("recombinationPlotFromVCF: putting parents first\n")
    }
    GTsInt <- cbind(
      GTsInt[, parentalIDs],
      GTsInt[, -which(dimnames(GTsInt)[[2]] %in% parentalIDs)]
    )
  }
  if(removeMissingInParents[1] == "either") {
    if(verbose) {
      cat("recombinationPlotFromVCF: removing unwanted variants\n")
    }
    if(shouldRemoveMendelianErrors) {
      rowsToUse <- (!GTsInt[, 1]==0 & !GTsInt[, 2]==0) & !(GTsInt[, 1] == GTsInt[, 2])
    } else {
      rowsToUse <- !GTsInt[, 1]==0 & !GTsInt[, 2]==0
    }
    GTsInt <- GTsInt[rowsToUse, ]
  }
  if(removeMissingInParents[1] == "both") {
    if(verbose) {
      cat("recombinationPlotFromVCF: removing unwanted variants\n")
    }
    if(shouldRemoveMendelianErrors) {
      rowsToUse <- (!GTsInt[, 1]==0 | !GTsInt[, 2]==0) & !(GTsInt[, 1] == GTsInt[, 2])
    } else {
      rowsToUse <- !GTsInt[, 1]==0 | !GTsInt[, 2]==0
    }
    GTsInt <- GTsInt[rowsToUse, ]
  }

  if(verbose) {
    cat("recombinationPlotFromVCF: creating plot\n")
  }
  recombinationPlot(
    convertGTsIntToParentBasedGTs(
      GTsInt
    ),
    ...
  )
}
