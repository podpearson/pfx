# recombinationPlotFromVCF.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################

# pdf("MAL4_7g8xGb4_release0.1.pdf", width=14, height=4)
# recombinationPlotFromVCF(
#   "/data/galton/users/rpearson/crossesTesting/release/7g8xGb4-qcPlusSamples-0.1.vcf.gz",
#   "MAL4"
# )
# dev.off()
#
# pdf("MAL4_7g8xGb4_Jiang.pdf", width=14, height=4)
# recombinationPlotFromVCF(
#   "/data/malariagen2/users/rpearson/pfCrosses/externalData/gb-2011-12-4-r33-s3.vcf.gz",
#   "MAL4",
#   GTsToIntMapping=c("7"=1, "G"=2)
# )
# dev.off()
# 
# pdf("MAL4_7g8xGb4_Zam_SNPs_20_200.pdf", width=14, height=4)
# recombinationPlotFromVCF(
#   "/data/galton/users/rpearson/zam/delivery/plasmodium/7g8_gb4_wk_flow_I_combined_BC_calls_at_all_k.decomp.vcf",
#   "MAL4",
#   GTsToIntMapping=c("0/0"=1, "1/1"=2),
#   keepPASSvariantsOnly=TRUE,
#   additionalInfoFilters     = list("SVTYPE" = list(operator="%in%", value="SNP")),
#   additionalGenotypeFilters     = list("GT_CONF" = list(operator="<=", value=20), "SITE_CONF" = list(operator="<=", value=200))
# )
# dev.off()

recombinationPlotFromVCF <- function(
  vcfFilename                 = "/data/galton/users/rpearson/crossesTesting/release/7g8xGb4-qcPlusSamples-0.1.vcf.gz",
  chromosome                  = "MAL4",
  parentalIDs                 = NULL, # this default will use IDs of first two samples in VCF. This can be over-ridden by supplying IDs, e.g. recombinationPlotFromVCF(parentalIDs=c("ERR027099", "ERR027100"))
  shouldRemoveInvariant       = TRUE,
  removeMissingInParents      = c("either", "both", "ignore"),
  shouldRemoveMendelianErrors = TRUE,
  regionsMask                 = varRegions_v2(), # will remove any variants in these regions. Set to NULL if you don't want to mask any variants out in this way
  keepPASSvariantsOnly        = FALSE,
  filtersToRemove             = c("Heterozygous", "Missingness"),     # set to NULL if you want to keep all variants. Also, this is overridden if keepPASSvariantsOnly==TRUE
  samplesToRemove             = NULL,
  additionalInfoFilters       = NULL,
#  additionalInfoFilters     = list(
#    "SVTYPE" = list(operator="%in%", value="SNP")
#  ),
  additionalGenotypeFilters   = NULL,
#  additionalGenotypeFilters     = list(
#    "GT_CONF" = list(operator="<=", value=20),
#    "SITE_CONF" = list(operator="<=", value=200)
#  )
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0, "./."=0, "2"=0, "3"=0), # "./." is needed as sometimes this is output by GATK's UG (presumably a bug). "2", "3", needed for the case of multi-allelic sites
  parametersList              = list(
    "7g8xGb4" = list(
      vcfFilename                 = vcfFilename,
      parentalIDs                 = parentalIDs,
      shouldRemoveInvariant       = shouldRemoveInvariant,
      removeMissingInParents      = removeMissingInParents,
      shouldRemoveMendelianErrors = shouldRemoveMendelianErrors,
      regionsMask                 = regionsMask,
      keepPASSvariantsOnly        = keepPASSvariantsOnly,
      filtersToRemove             = filtersToRemove,
      samplesToRemove             = samplesToRemove,
      additionalInfoFilters       = additionalInfoFilters,
      additionalGenotypeFilters   = additionalGenotypeFilters,
      GTsToIntMapping             = c("0"=1, "1"=2, "."=0, "./."=0, "2"=0, "3"=0) # "./." is needed as sometimes this is output by GATK's UG (presumably a bug). "2", "3", needed for the case of multi-allelic sites
    )
  ),
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
      keepPASSvariantsOnly        = keepPASSvariantsOnly,
      filtersToRemove             = filtersToRemove,
      samplesToRemove             = samplesToRemove,
      additionalInfoFilters       = additionalInfoFilters,
      additionalGenotypeFilters   = additionalGenotypeFilters,
      overwriteExisting           = TRUE,
      saveAsRobjectFile           = FALSE
    ),
    GTsToIntMapping             = GTsToIntMapping
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
