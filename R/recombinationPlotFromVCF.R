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
#   GTsToIntMapping=c("7"=1, "G"=2),
#   sampleIDmappings=createSampleIDmappings(
#     "/data/malariagen2/users/rpearson/pfCrosses/externalData/gb-2011-12-4-r33-s3.vcf.gz",
#     shouldUseSampleAnnotation = FALSE,
#     sampleDuplicates        = c(
#       "_7G8" = "7G8_JH6",
#       "JH6"  = "7G8_JH6",
#       "JF6"  = "JF6_KC2",
#       "KC2"  = "JF6_KC2",
#       "AUD"  = "AUD_LC12",
#       "LC12" = "AUD_LC12",
#       "KB8"  = "KB8_KC5_NH11",
#       "KC5"  = "KB8_KC5_NH11",
#       "NH11" = "KB8_KC5_NH11",
#       "D2"   = "D2_TF1",
#       "TF1"  = "D2_TF1"
#     )
#   )
# )
# dev.off()
# 
# pdf("MAL4_7g8xGb4_Zam.pdf", width=14, height=4)
# recombinationPlotFromVCF(
#   "/data/galton/users/rpearson/zam/delivery/plasmodium/7g8_gb4_wk_flow_I_combined_BC_calls_at_all_k.decomp.vcf.gz",
#   "MAL4",
#   GTsToIntMapping=c("0/0"=1, "1/1"=2, "."=0),
#   keepPASSvariantsOnly=TRUE
# )
# dev.off()
# pdf("MAL4_7g8xGb4_Zam_SNPs_20_200.pdf", width=14, height=4)
# recombinationPlotFromVCF(
#   "/data/galton/users/rpearson/zam/delivery/plasmodium/7g8_gb4_wk_flow_I_combined_BC_calls_at_all_k.decomp.vcf.gz",
#   "MAL4",
#   GTsToIntMapping=c("0/0"=1, "1/1"=2, "."=0),
#   keepPASSvariantsOnly=TRUE,
#   additionalInfoFilters     = list("SVTYPE" = list(operator="%in%", value="SNP")),
#   additionalGenotypeFilters     = list("GT_CONF" = list(operator="<=", value=20), "SITE_CONF" = list(operator="<=", value=200))
# )
# dev.off()
# pdf("MAL4_7g8xGb4_Zam_SNPs_50_1000.pdf", width=14, height=4)
# GTsInt <- recombinationPlotFromVCF(
#   "/data/galton/users/rpearson/zam/delivery/plasmodium/7g8_gb4_wk_flow_I_combined_BC_calls_at_all_k.decomp.vcf.gz",
#   "MAL4",
#   GTsToIntMapping=c("0/0"=1, "1/1"=2, "."=0),
#   keepPASSvariantsOnly=TRUE,
#   additionalInfoFilters     = list("SVTYPE" = list(operator="%in%", value="SNP")),
#   additionalGenotypeFilters     = list("GT_CONF" = list(operator="<=", value=50), "SITE_CONF" = list(operator="<=", value=1000))
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
  genoToLoad                  = c("GT", "AD", names(additionalGenotypeFilters)),
  GTsToIntMapping             = c("0"=1, "0/0"=1, "0|0"=1, "1"=2, "1/1"=2, "1|1"=2, "."=0, "./."=0, "./."=0, "2"=0, "3"=0, "0/1"=0, "1/0"=0, "0|1"=0, "1|0"=0), # "./." is needed as sometimes this is output by GATK's UG (presumably a bug). "2", "3", needed for the case of multi-allelic sites
#  GTsToIntMapping             = c("0"=1, "1"=2, "."=0, "./."=0, "2"=0, "3"=0), # "./." is needed as sometimes this is output by GATK's UG (presumably a bug). "2", "3", needed for the case of multi-allelic sites
  sampleIDmappings            = NULL,
#  sampleIDmappings            = createSampleIDmappings(vcfFilename),
  linePositions               = c(2),
  parametersList              = list(
    "7g8xGb4" = list(
      vcfFilename                 = vcfFilename,
      parentalIDs                 = parentalIDs,
      genoToLoad                  = genoToLoad,
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
  GTsIntList <- sapply(
    parametersList,
    function(parameters) {
      GTsInt <- genotypeCallsFromGTas012(
        createSingleChromosomeVariantSitesRdaFile(
          vcfFilename                 = parameters[["vcfFilename"]],
          chromosome                  = chromosome,
          genoToLoad                  = parameters[["genoToLoad"]],
          shouldRemoveInvariant       = parameters[["shouldRemoveInvariant"]],
          regionsMask                 = parameters[["regionsMask"]],
          keepPASSvariantsOnly        = parameters[["keepPASSvariantsOnly"]],
          filtersToRemove             = parameters[["filtersToRemove"]],
          samplesToRemove             = parameters[["samplesToRemove"]],
          additionalInfoFilters       = parameters[["additionalInfoFilters"]],
          additionalGenotypeFilters   = parameters[["additionalGenotypeFilters"]],
          shouldAnnotateUsingExternal = FALSE,
          overwriteExisting           = TRUE,
          saveAsRobjectFile           = FALSE
        ),
        GTsToIntMapping             = GTsToIntMapping
      )
#      browser()
      if(is.null(parentalIDs)) {
        parentalIDs <<- dimnames(GTsInt)[[2]][1:2]
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
      if(!is.null(sampleIDmappings)) {
        return(reorderSamples(GTsInt, parentalIDs, sampleIDmappings))
      } else {
        return(list(GTsInt=GTsInt, linePositions=linePositions, parentalIDs=parentalIDs))
      }
    },
    USE.NAMES=TRUE,
    simplify=FALSE
  )
#  GTsInt <- genotypeCallsFromGTas012(
#    createSingleChromosomeVariantSitesRdaFile(
#      vcfFilename                 = vcfFilename,
#      chromosome                  = chromosome,
#      genoToLoad                  = genoToLoad,
#      shouldRemoveInvariant       = shouldRemoveInvariant,
#      regionsMask                 = regionsMask,
#      keepPASSvariantsOnly        = keepPASSvariantsOnly,
#      filtersToRemove             = filtersToRemove,
#      samplesToRemove             = samplesToRemove,
#      additionalInfoFilters       = additionalInfoFilters,
#      additionalGenotypeFilters   = additionalGenotypeFilters,
#      overwriteExisting           = TRUE,
#      saveAsRobjectFile           = FALSE
#    ),
#    GTsToIntMapping             = GTsToIntMapping
#  )
#  if(is.null(parentalIDs)) {
#    parentalIDs <- dimnames(GTsInt)[[2]][1:2]
#  }
#  if(!identical(parentalIDs, dimnames(GTsInt)[[2]][1:2])) {
#    if(verbose) {
#      cat("recombinationPlotFromVCF: putting parents first\n")
#    }
#    GTsInt <- cbind(
#      GTsInt[, parentalIDs],
#      GTsInt[, -which(dimnames(GTsInt)[[2]] %in% parentalIDs)]
#    )
#  }
#  if(removeMissingInParents[1] == "either") {
#    if(verbose) {
#      cat("recombinationPlotFromVCF: removing unwanted variants\n")
#    }
#    if(shouldRemoveMendelianErrors) {
#      rowsToUse <- (!GTsInt[, 1]==0 & !GTsInt[, 2]==0) & !(GTsInt[, 1] == GTsInt[, 2])
#    } else {
#      rowsToUse <- !GTsInt[, 1]==0 & !GTsInt[, 2]==0
#    }
#    GTsInt <- GTsInt[rowsToUse, ]
#  }
#  if(removeMissingInParents[1] == "both") {
#    if(verbose) {
#      cat("recombinationPlotFromVCF: removing unwanted variants\n")
#    }
#    if(shouldRemoveMendelianErrors) {
#      rowsToUse <- (!GTsInt[, 1]==0 | !GTsInt[, 2]==0) & !(GTsInt[, 1] == GTsInt[, 2])
#    } else {
#      rowsToUse <- !GTsInt[, 1]==0 | !GTsInt[, 2]==0
#    }
#    GTsInt <- GTsInt[rowsToUse, ]
#  }

  if(verbose) {
    cat("recombinationPlotFromVCF: creating plot\n")
  }
#  browser()
  sapply(
    GTsIntList,
    function(GTsInt) {
      recombinationPlot(
        convertGTsIntToParentBasedGTs(
          GTsInt[["GTsInt"]],
          IDparent1 = GTsInt[["parentalIDs"]][1],
          IDparent2 = GTsInt[["parentalIDs"]][2]
        ),
        linePositions = GTsInt[["linePositions"]],
        ...
      )
    }
  )
  return(GTsIntList)
}
