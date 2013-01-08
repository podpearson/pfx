# filterEvaluationSeries.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


filterEvaluationSeries <- function(
  cross                       = "3d7_hb3",
  variantType                 = "snps",
  analysisDirectory           = "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20120107/gatk",
  filters = list(
    "homopolymer5Proximity" = list(column="homopolymer5Proximity", operator="%in%", value=1),
    "homopolymer15Proximity" = list(column="homopolymer15Proximity", operator="%in%", value=0:20),
    "MQ0" = list(column="MQ0", operator=">", value=0),
    "TRF" = list(column="RepeatPeriod1", operator=">", value=0),
    "ReadPosRankSum" = list(column="ReadPosRankSum", operator="<", value=-2, filterOutNAs=TRUE),
    "scaledDepthSD" = list(column="scaledDepthSD", operator=">", value=0.5),
    "SoftClipped" = list(column="SoftClipped", operator=">", value=0.1),
    "UQ35" = list(column="UQ", operator=">", value=35),
    "QUAL12000" = list(column="QUAL", operator="<", value=12000)
  )
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
  load(file.path(analysisDirectory, cross, variantType, paste(cross, ".vcfCoreFinalSamples.rda", sep="")))
  
  filterResults <- sapply(
    seq(along=filters),
    function(numberOfFilters) {
      filtersToUse <- filters[1:numberOfFilters]
      currentFiltersResults <- evaluateFilters(
        vcfCoreFinalSamples,
        plotFilestem=file.path(analysisDirectory, cross, variantType, "evaluateFilters"),
        additionalInfoFilters = filtersToUse,
        sampleDuplicates=initialSampleQCresults[["sampleDuplicates"]]
      )
      return(currentFiltersResults)
    },
    simplify=FALSE
  )
  filterResultsDF <- do.call(rbind, filterResults)
  filterResultsDFfilename <- file.path(analysisDirectory, cross, variantType, paste(cross, ".filterResultsDF.rda", sep=""))
  save(filterResultsDF, file=filterResultsDFfilename)
  return(filterResultsDF)
}
