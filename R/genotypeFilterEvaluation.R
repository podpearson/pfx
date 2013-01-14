# genotypeFilterEvaluation.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


#  genotypeFilterEvaluation2_3d7_hb3_snps <- genotypeFilterEvaluation("3d7_hb3", "snps", monomorphicSkipChromosomes  = "Pf3D7_13_v3")
#  genotypeFilterEvaluation2_7g8_gb4_snps <- genotypeFilterEvaluation("7g8_gb4", "snps")
#  genotypeFilterEvaluation2_hb3_dd2_snps <- genotypeFilterEvaluation("hb3_dd2", "snps")
#  genotypeFilterEvaluation2_3d7_hb3_indels <- genotypeFilterEvaluation("3d7_hb3", "indels", monomorphicSkipChromosomes  = "Pf3D7_13_v3", minMeanMAFtoConsiderContam=0.02)
#  genotypeFilterEvaluation2_7g8_gb4_indels <- genotypeFilterEvaluation("7g8_gb4", "indels", minMeanMAFtoConsiderContam=0.02)
#  genotypeFilterEvaluation2_hb3_dd2_indels <- genotypeFilterEvaluation("hb3_dd2", "indels", minMeanMAFtoConsiderContam=0.02)

genotypeFilterEvaluation <- function(
  cross                       = "3d7_hb3",
  variantType                 = "snps",
  analysisDirectory           = "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20120107/gatk",
#  filters                     = NULL,
  filters=list(
#      "biallelic" = list(column="biallelic", operator="==", value=TRUE),
    "meanMAF0.2" = list(column="meanMAF", operator=">", value=0.2, filterOutNAs=TRUE),
    "QD" = list(column="QUALbyDP", operator=">", value=33),
    "missingness1" = list(column="missingness", operator=">", value=1),
    "numFilteredGenotypes2" = list(column="numFilteredGenotypes", operator=">", value=2),
    "MQ0Fraction" = list(column="MQ0Fraction", operator=">", value=0.01),
    "MateOtherChrom" = list(column="MateOtherChrom", operator=">", value=0.1),
    "SoftClipped" = list(column="MateOtherChrom", operator=">", value=0.2),
    "ProperPair" = list(column="ProperPair", operator="<", value=0.9),
    "MQ" = list(column="MQ", operator="<", value=38),
    "HRun" = list(column="HRun", operator=">", value=4, filterOutNAs=TRUE),
    "FS" = list(column="FS", operator=">", value=600)
  ),
#  regionsMask                 = NULL,
  regionsMask                 = varRegions_v3(),
  setMonomorphicProgenyFilter = FALSE,
  monomorphicSkipChromosomes  = NULL,
  GQthresholds                = c(99, 50, 5),
  DPthresholds                = c(10, 5, 1),
  MAFthresholds               = c(0.1, 0, 0.02, 0.05, 0.2),
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
  minMeanMAFtoConsiderContam  = 0.01,
  shouldReturnVcfOnly         = FALSE,
  shouldUseExistingRda        = FALSE
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
  vcfAnnotatedBestReplicateSamplesFilename <- file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfAnnotatedBestReplicateSamples.rda", sep=""))
  vcfAnnotatedUncontaminatedSamplesFilename <- file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfAnnotatedUncontaminatedSamples.rda", sep=""))
  if(shouldUseExistingRda && file.exists(vcfAnnotatedFinalSamplesFilename) && file.exists(vcfAnnotatedBestReplicateSamplesFilename) && file.exists(vcfAnnotatedUncontaminatedSamplesFilename)) {
    load(vcfAnnotatedFinalSamplesFilename)
    load(vcfAnnotatedBestReplicateSamplesFilename)
    load(vcfAnnotatedUncontaminatedSamplesFilename)
  } else {
    if(is.null(regionsMask)) {
      load(file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfVariant.rda", sep="")))
    } else {
      load(file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfCoreFinalSamples.rda", sep="")))
      vcfVariant <- vcfCoreFinalSamples
      rm(vcfCoreFinalSamples)
      gc()
    }
    finalSamples <- setdiff(dimnames(vcfVariant)[[2]], initialSampleQCresults[["qcFailedSamples"]])
    vcfAnnotatedFinalSamples <-  annotateVcf(vcfVariant[, finalSamples])
    save(vcfAnnotatedFinalSamples, file=vcfAnnotatedFinalSamplesFilename)
#    bestReplicateSamples <- setdiff(dimnames(vcfVariant)[[2]], initialSampleQCresults[["uniqueSamples"]])
    bestReplicateSamples <- setdiff(initialSampleQCresults[["uniqueSamples"]], initialSampleQCresults[["qcFailedSamples"]])
    vcfAnnotatedBestReplicateSamples <-  annotateVcf(vcfVariant[, bestReplicateSamples])
    save(vcfAnnotatedBestReplicateSamples, file=vcfAnnotatedBestReplicateSamplesFilename)
    if(!exists("vcfInitialFiltered")) {
      load(file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfInitialFiltered.rda", sep="")))
    }
    vcfInitialFilteredPASS <- vcfInitialFiltered[filt(vcfInitialFiltered)=="PASS"]
    RefReads <- matrix(
      sapply(geno(vcfInitialFilteredPASS)[["AD"]], function(x) x[1]),
      ncol=dim(geno(vcfInitialFilteredPASS)[["AD"]])[2],
      dimnames=dimnames(geno(vcfInitialFilteredPASS)[["AD"]])
    )
    FirstAltReads <- matrix(
      sapply(geno(vcfInitialFilteredPASS)[["AD"]], function(x) x[2]),
      ncol=dim(geno(vcfInitialFilteredPASS)[["AD"]])[2],
      dimnames=dimnames(geno(vcfInitialFilteredPASS)[["AD"]])
    )
    MAF <- pmin(RefReads, FirstAltReads)/(RefReads+FirstAltReads)
    meanMAFperSample <- colMeans(MAF, na.rm = TRUE)
    uncontaminatedSamples <- which(meanMAFperSample < minMeanMAFtoConsiderContam)
    vcfAnnotatedUncontaminatedSamples <-  annotateVcf(vcfVariant[, uncontaminatedSamples])
    save(vcfAnnotatedUncontaminatedSamples, file=vcfAnnotatedUncontaminatedSamplesFilename)
  }
  
#  browser()
#  
#  filterResults <- evaluateGenotypeFilters(
#    vcfAnnotatedFinalSamples,
##    vcfCoreFinalSamples,
#    plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", sep=".")),
#    additionalInfoFilters       = filters,
#    regionsMask                 = regionsMask,
#    genotypeFiltersList         = genotypeFiltersList,
#    setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
#    monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
#    maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
#    sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
#    shouldReturnVcfOnly         = shouldReturnVcfOnly
#  )
  filterResultsFinalSamples <- evaluateGenotypeFilters(
    vcfAnnotatedFinalSamples,
#    vcfCoreFinalSamples,
    plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", "final", sep=".")),
    additionalInfoFilters       = filters,
    regionsMask                 = regionsMask,
    genotypeFiltersList         = genotypeFiltersList,
    setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
    monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
    maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
    sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
    shouldReturnVcfOnly         = shouldReturnVcfOnly
  )
  names(filterResultsFinalSamples) <- paste("FinalSamples", names(filterResultsFinalSamples))
  filterResultsBestReplicateSamples <- evaluateGenotypeFilters(
    vcfAnnotatedBestReplicateSamples,
#    vcfCoreFinalSamples,
    plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", "bestReplicate", sep=".")),
    additionalInfoFilters       = filters,
    regionsMask                 = regionsMask,
    genotypeFiltersList         = genotypeFiltersList,
    setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
    monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
    maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
    sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
    shouldReturnVcfOnly         = shouldReturnVcfOnly
  )
  names(filterResultsBestReplicateSamples) <- paste("BestReplicate", names(filterResultsBestReplicateSamples))
  filterResultsUncontaminatedSamples <- evaluateGenotypeFilters(
    vcfAnnotatedUncontaminatedSamples,
#    vcfCoreFinalSamples,
    plotFilestem                = file.path(analysisDirectory, cross, variantType, paste(cross, variantType, "evaluateGenotypeFilters", "uncontaminated", sep=".")),
    additionalInfoFilters       = filters,
    regionsMask                 = regionsMask,
    genotypeFiltersList         = genotypeFiltersList,
    setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
    monomorphicSkipChromosomes  = monomorphicSkipChromosomes,
    maxNumFilteredGenotypes     = maxNumFilteredGenotypes,
    sampleDuplicates            = initialSampleQCresults[["sampleDuplicates"]],
    shouldReturnVcfOnly         = shouldReturnVcfOnly
  )
  names(filterResultsUncontaminatedSamples) <- paste("UncontaminatedSamples", names(filterResultsUncontaminatedSamples))
  filterResults <- cbind(filterResultsFinalSamples, filterResultsBestReplicateSamples, filterResultsUncontaminatedSamples)
  return(filterResults)
}
