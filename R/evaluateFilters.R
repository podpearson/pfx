# evaluateFilters.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


evaluateFilters <- function(
  vcf,
  additionalInfoFilters = list(
    "LowQD" = list(column="QD", operator="<=", value=36),
    "HighSB" = list(column="SB", operator=">=", value=-6000)
  ),
  regionsMask                 = varRegions_v3(),
  regionsMaskFilterName       = "InVarRegion",
  plotFilestem                = "evaluateFilters",
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0, "./."=0),
  sampleIDcolumn              = "ena_run_accession",
  sampleIDmappingsColumn      = sampleIDcolumn,
  sampleDuplicates            = NULL,
  shouldRecalculateDepthSD    = TRUE
) {
  vcfFiltered <- filterVcf(
    setVcfFilters(
      vcf,
      regionsMask                 = regionsMask,
      additionalInfoFilters       = additionalInfoFilters
    ),
    filtersToRemove = c(regionsMaskFilterName, names(additionalInfoFilters))
  )
  if(shouldRecalculateDepthSD) {
    values(info(vcfFiltered))[["scaledDepthSD"]] <- calcuateScaledDepthSD(vcfFiltered)
  }
  qcFilteringPlots(vcfFiltered, plotFilestem=paste(c(plotFilestem, regionsMaskFilterName, names(additionalInfoFilters)), collapse="."), shouldCreateErrorRateBySites=FALSE)
  mgRecombinations <- recombinationPoints(vcfFiltered, shouldCharacterise=FALSE, GTsToIntMapping=GTsToIntMapping)
  recombinationsPerSample <- rev(
    sapply(
      mgRecombinations[["sampleLevelResults"]],
      function(x) {
        sum(sapply(x, length))
      }
    )
  )
  recombinationPlotSeries(
    vcfFiltered,
    plotFilestem                = paste(c(plotFilestem, regionsMaskFilterName, names(additionalInfoFilters)), collapse="."),
    filters                     = NULL,
    sampleIDcolumn              = sampleIDcolumn,
    sampleIDmappingsColumn      = sampleIDmappingsColumn,
    sampleDuplicates            = sampleDuplicates
  )
#  REFs <- as.character(ref(vcfCoreFinalSamples))
#  ALTs <- as.character(unlist(alt(vcfCoreFinalSamples)))[do.call(c, sapply(elementLengths(alt(vcfCoreFinalSamples)), seq))==1]
#  transitions <- REFs=="A" & ALTs=="G" | REFs=="G" & ALTs=="A" | REFs=="C" & ALTs=="T" | REFs=="T" & ALTs=="C" 
#  transversions <- REFs=="A" & ALTs %in% c("T", "C") | REFs=="G" & ALTs %in% c("T", "C") | REFs=="C" & ALTs %in% c("A", "G") | REFs=="T" & ALTs %in% c("A", "G")
#  titvRatio <- length(which(transitions)) / length(which(transversions))
  titvRatio <- titv(vcfFiltered)
  titvRatioExcludingAT <- titv(vcfFiltered, FALSE)

  returnDF <- data.frame(
    totalRecombinations = sum(recombinationsPerSample),
    medianRecombinationsPerSample = median(recombinationsPerSample),
    numberOfVariants = dim(vcfFiltered)[1],
    numberOfMendelianErrors = length(which(values(info(vcfFiltered))[["MendelianErrors"]] > 0)),
    numberOfSegregatingSites = length(which(values(info(vcfFiltered))[["SEGREGATING"]])),
    titvRatio = titvRatio,
    titvRatioExcludingAT = titvRatioExcludingAT,
    row.names = paste(c(regionsMaskFilterName, names(additionalInfoFilters)), collapse=".")
  )
  save(returnDF, file=paste(paste(c(plotFilestem, regionsMaskFilterName, names(additionalInfoFilters)), collapse="."), "returnDF.rda", sep="."))
  
  return(returnDF)
}
