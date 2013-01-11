# genotypeFilterEvaluation.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################



genotypeFilterEvaluation <- function(
  cross                       = "3d7_hb3",
  variantType                 = "snps",
  analysisDirectory           = "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20120107/gatk",
  filters                     = NULL,
  regionsMask                 = NULL,
#  regionsMask                 = varRegions_v3(),
  setMonomorphicProgenyFilter = TRUE,
  GQthresholds                = c(99, 50, 5),
  DPthresholds                = c(10, 5, 1),
  MAFthresholds               = c(0, 0.02, 0.05, 0.1, 0.2),
#  GQthresholds                = c(99, seq(95, 0, -5)),
#  DPthresholds                = seq(20, 1, -1),
#  MAFthresholds               = c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5),
  GQthresholdDefault          = 99,
  DPthresholdDefault          = 10,
  MAFthresholdDefault         = 0.1,
  genotypeFiltersList         = c(
    lapply(
      MAFthresholds,
      function(MAFthreshold) {
        list(
          "LowGQ" = list(column="GQ", operator="<", value=GQthresholdDefault, filterOutNAs=TRUE),
          "LowDP" = list(column="DP", operator="<", value=DPthresholdDefault, filterOutNAs=TRUE),
          "HighMAF" = list(column="MAF", operator=">", value=MAFthreshold, filterOutNAs=TRUE)
        )
      }
    ),
    lapply(
      DPthresholds,
      function(DPthreshold) {
        list(
          "LowGQ" = list(column="GQ", operator="<", value=GQthresholdDefault, filterOutNAs=TRUE),
          "LowDP" = list(column="DP", operator="<", value=DPthreshold, filterOutNAs=TRUE),
          "HighMAF" = list(column="MAF", operator=">", value=MAFthresholdDefault, filterOutNAs=TRUE)
        )
      }
    ),
    lapply(
      GQthresholds,
      function(GQthreshold) {
        list(
          "LowGQ" = list(column="GQ", operator="<", value=GQthreshold, filterOutNAs=TRUE),
          "LowDP" = list(column="DP", operator="<", value=DPthresholdDefault, filterOutNAs=TRUE),
          "HighMAF" = list(column="MAF", operator=">", value=MAFthresholdDefault, filterOutNAs=TRUE)
        )
      }
    )
  ),
  maxNumFilteredGenotypes     = 2,
  shouldReturnVcfOnly         = FALSE
)
{
  initialSampleQCresultsFilename <- file.path(analysisDirectory, cross, variantType, paste(cross, ".initialSampleQCresults.rda", sep=""))
  if(file.exists(initialSampleQCresultsFilename)) {
    load(initialSampleQCresultsFilename)
  } else {
    load(file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfInitialFiltered.rda", sep="")))
    vcfSegregating <- filterVcf(vcfInitialFiltered, keepPASSvariantsOnly=TRUE)
    initialSampleQCresults <- sampleQC(
      vcfSegregating,
      discordanceThreshold=1000,
      shouldCreatePlots=FALSE,
      sampleIDcolumn="ena_run_accession",
      sampleIDmappingsColumn="ena_run_accession"
    )
    save(initialSampleQCresults, file=initialSampleQCresultsFilename)
  }
  vcfAnnotatedFinalSamplesFilename <- file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfAnnotatedFinalSamples.rda", sep=""))
  if(file.exists(vcfAnnotatedFinalSamplesFilename)) {
    load(vcfAnnotatedFinalSamplesFilename)
  } else {
    load(file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfVariant.rda", sep="")))
    finalSamples <- setdiff(dimnames(vcfVariant)[[2]], initialSampleQCresults[["qcFailedSamples"]])
    vcfAnnotatedFinalSamples <-  annotateVcf(vcfVariant[, finalSamples])
    save(vcfAnnotatedFinalSamples, file=vcfAnnotatedFinalSamplesFilename)
  }
  
  filterResults <- evaluateGenotypeFilters(
    vcfAnnotatedFinalSamples,
#    vcfCoreFinalSamples,
    plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", sep=".")),
    additionalInfoFilters       = filters,
    regionsMask                 = regionsMask,
    genotypeFiltersList         = genotypeFiltersList,
    setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
    maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
    sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
    shouldReturnVcfOnly         = shouldReturnVcfOnly
  )
  return(filterResults)
}
