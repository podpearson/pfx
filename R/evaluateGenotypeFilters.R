# evaluateGenotypeFilters.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


evaluateGenotypeFilters <- function(
  vcf,
  additionalInfoFilters = NULL,
#  additionalInfoFilters = list(
#    "LowQD" = list(column="QD", operator="<=", value=36),
#    "HighSB" = list(column="SB", operator=">=", value=-6000)
#  ),
  regionsMask                 = NULL,
  regionsMaskFilterName       = "noMask",
#  regionsMask                 = varRegions_v3(),
#  regionsMaskFilterName       = "InVarRegion",
  shouldSetMultiallelicFilter = TRUE,
  shouldSetNonSegregatingFilt = FALSE,
  shouldSetMissingInParentFilt= TRUE,
  setMonomorphicProgenyFilter = TRUE,
  monomorphicSkipChromosomes  = NULL,
  plotFilestem                = "evaluateGenotypeFilters",
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0, "./."=0),
  possibleMissingValues       = c(".", "./.", ".|."),
  parentalIDs                 = dimnames(vcf)[[2]][1:2],
  parentalIDindexes           = which(dimnames(vcf)[[2]] %in% parentalIDs),
  IDsToRemoveFromDuplicates   = "ERR027105",
  sampleIDcolumn              = "ena_run_accession",
  sampleIDmappingsColumn      = sampleIDcolumn,
  sampleDuplicates            = NULL,
  shouldRecalculateDepthSD    = FALSE,
  shouldCalculateExtraQUAL    = FALSE,
  shouldFilterGenotypes       = TRUE,
  shouldCreateQCFilteringPlots= FALSE,
  shouldCreateRecombPlots     = TRUE,
#  genotypeFilters             = list(
#    "LowGQ" = list(column="GQ", operator="<", value=99, filterOutNAs=TRUE),
#    "LowDP" = list(column="DP", operator="<", value=10, filterOutNAs=TRUE),
#    "HighMAF" = list(column="MAF", operator=">", value=0.1, filterOutNAs=TRUE)
#  ),
  GQthresholds                = c(99, 50, 5),
  DPthresholds                = seq(10, 5, 1),
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
  errorVariables              = c("MendelianErrors", "numSingleSNPhaplotypes"),
  shouldSetHaplotypeLengths   = TRUE,
  shouldReturnVcfOnly         = FALSE
) {
  resultsList <- lapply(
    genotypeFiltersList,
    function(genotypeFilters) {
      cat("evaluateGenotypeFilters:", genotypeFilters[["LowGQ"]][["value"]], genotypeFilters[["LowDP"]][["value"]], genotypeFilters[["HighMAF"]][["value"]], "\n")
      theseFiltersPlotFilestem <- paste(
        c(
          plotFilestem,
          genotypeFilters[["LowGQ"]][["value"]],
          genotypeFilters[["LowDP"]][["value"]],
          genotypeFilters[["HighMAF"]][["value"]],
          sapply(additionalInfoFilters, function(x) x[["value"]])
        ),
        collapse="."
      )
      if(shouldFilterGenotypes) {
        vcfFiltered <- filterGenotypes(
          vcf,
          genotypeFilters             = genotypeFilters
    #      shouldSetNonSegregatingFilt = shouldSetNonSegregatingFilt,
    #      regionsMask                 = regionsMask,
    #      maxNumFilteredGenotypes     = maxNumFilteredGenotypes
        )
      }
      vcfFiltered <- filterVcf(
        setVcfFilters(
          vcfFiltered,
          regionsMask                 = regionsMask,
          additionalInfoFilters       = additionalInfoFilters,
          shouldSetMultiallelicFilter = shouldSetMultiallelicFilter,
          shouldSetNonSegregatingFilt = shouldSetNonSegregatingFilt,
          shouldSetMissingInParentFilt=shouldSetMissingInParentFilt,
          setMonomorphicProgenyFilter = setMonomorphicProgenyFilter,
          monomorphicSkipChromosomes  = monomorphicSkipChromosomes
        ),
        keepPASSvariantsOnly = TRUE
      )
      browser()
      filterColumns <- sapply(additionalInfoFilters, function(x) x[["column"]])
      if(shouldRecalculateDepthSD) {
        currentInfo <- values(info(vcfFiltered))
        currentInfo[["scaledDepthSD"]] <- calcuateScaledDepthSD(vcfFiltered)
        info(vcfFiltered) <- currentInfo
        if("scaledDepthSD" %in% filterColumns) {
          vcfFiltered <- filterVcf(
            setVcfFilters(
              vcfFiltered,
              regionsMask                 = regionsMask,
              additionalInfoFilters       = additionalInfoFilters
            ),
            filtersToRemove = c(regionsMaskFilterName, names(additionalInfoFilters))
          )
        }
      }
      if(shouldCalculateExtraQUAL) {
        currentInfo <- values(info(vcfFiltered))
        currentInfo[["QUALbyDP"]] = qual(vcfFiltered)/values(info(vcfFiltered))[["DP"]]
        currentInfo[["QUALperSample"]] = qual(vcfFiltered)/dim(vcfFiltered)[2]
        info(vcfFiltered) <- currentInfo
        if("QUALbyDP" %in% filterColumns || "QUALperSample" %in% filterColumns) {
          vcfFiltered <- filterVcf(
            setVcfFilters(
              vcfFiltered,
              regionsMask                 = regionsMask,
              additionalInfoFilters       = additionalInfoFilters
            ),
            filtersToRemove = c(regionsMaskFilterName, names(additionalInfoFilters))
          )
        }
      }
      if(shouldSetHaplotypeLengths) {
        vcfFiltered <- setHaplotypeLengths(vcfFiltered)
      }
    #  browser()
      if(shouldReturnVcfOnly) {
        return(vcfFiltered)
      }
      if(shouldCreateQCFilteringPlots) {
        sapply(
          errorVariables,
          function(errorVariable) {
            qcFilteringPlots(
              vcfFiltered,
              plotFilestem=theseFiltersPlotFilestem,
      #        plotFilestem=paste(c(plotFilestem, regionsMaskFilterName, names(additionalInfoFilters), errorVariable), collapse="."),
              shouldCreateErrorRateBySites=FALSE,
              errorVariable=errorVariable
            )
          }
        )
      }
##      The following is buggy - not identify recombinations where there are missing genotypes
#      mgRecombinations <- recombinationPoints(vcfFiltered, shouldCharacterise=FALSE, GTsToIntMapping=GTsToIntMapping)
#      recombinationsPerSample <- rev(
#        sapply(
#          mgRecombinations[["sampleLevelResults"]],
#          function(x) {
#            sum(sapply(x, length))
#          }
#        )
#      )
      haplotypesPerSampleAndChromosomeMatrix <- haplotypesPerSample(vcfFiltered)
      haplotypesPerSample <- rowSums(haplotypesPerSampleAndChromosomeMatrix)
      singleSNPhaplotypesPerSampleAndChromosomeMatrix <- shortHaplotypesPerSample(vcfFiltered)
      singleSNPhaplotypesPerSample <- rowSums(singleSNPhaplotypesPerSampleAndChromosomeMatrix)
      putativeCrossoversPerSample <- haplotypesPerSample-14-(2*singleSNPhaplotypesPerSample)
      
      nonMissingGenotypesPerSample <- apply(
        geno(vcfFiltered)[["GT"]],
        2,
        function(x) {
          length(which(!(x %in% possibleMissingValues)))
        }
      )
      if(shouldCreateRecombPlots) {
        recombinationPlotSeries(
          vcfFiltered,
          plotFilestem                = theseFiltersPlotFilestem,
      #    plotFilestem                = paste(c(plotFilestem, regionsMaskFilterName, names(additionalInfoFilters)), collapse="."),
          filters                     = NULL,
          sampleIDcolumn              = sampleIDcolumn,
          sampleIDmappingsColumn      = sampleIDmappingsColumn,
          sampleDuplicates            = sampleDuplicates
        )
      }
    #  REFs <- as.character(ref(vcfCoreFinalSamples))
    #  ALTs <- as.character(unlist(alt(vcfCoreFinalSamples)))[do.call(c, sapply(elementLengths(alt(vcfCoreFinalSamples)), seq))==1]
    #  transitions <- REFs=="A" & ALTs=="G" | REFs=="G" & ALTs=="A" | REFs=="C" & ALTs=="T" | REFs=="T" & ALTs=="C" 
    #  transversions <- REFs=="A" & ALTs %in% c("T", "C") | REFs=="G" & ALTs %in% c("T", "C") | REFs=="C" & ALTs %in% c("A", "G") | REFs=="T" & ALTs %in% c("A", "G")
    #  titvRatio <- length(which(transitions)) / length(which(transversions))
      titvRatio <- titv(vcfFiltered)
      titvRatioExcludingAT <- titv(vcfFiltered, FALSE)
      
      sampleDuplicatesAsDF <- do.call(
        rbind,
        lapply(
          intersect(dimnames(vcfFiltered)[[2]], setdiff(unique(names(sampleDuplicates)), IDsToRemoveFromDuplicates)),
          function(sampleID) {
            duplicateIDs <- strsplit(sampleDuplicates[sampleID], "_")[[1]]
            duplicateIDs <- duplicateIDs[duplicateIDs != sampleID]
            duplicatePairs <- do.call(rbind, lapply(duplicateIDs, function(x) if(x>sampleID) data.frame(sampleID1=sampleID, sampleID2=x, stringsAsFactors=FALSE) else data.frame(sampleID1=character(0), sampleID2=character(0))))
          }
        )
      )
      duplicateDiscordanceMatrix <- duplicateDiscordances(vcfFiltered, sampleDuplicatesAsDF, possibleMissingValues)
      meanDuplicateDiscordanceRate <- mean(colSums(duplicateDiscordanceMatrix, na.rm=TRUE))
    
#      vcf_Pf3D7_01_v3 <- vcfFiltered[seqnames(vcfFiltered)=="Pf3D7_01_v3"]
#      seqlevels(vcf_Pf3D7_01_v3) <- "Pf3D7_01_v3"
#      mgRecombinations_Pf3D7_01_v3 <- recombinationPoints(vcf_Pf3D7_01_v3, shouldCharacterise=FALSE, GTsToIntMapping=GTsToIntMapping)
#      recombinationsPerSample_Pf3D7_01_v3 <- rev(
#        sapply(
#          mgRecombinations_Pf3D7_01_v3[["sampleLevelResults"]],
#          function(x) {
#            sum(sapply(x, length))
#          }
#        )
#      )
      returnDF <- data.frame(
        haplotypeParent1 = haplotypesPerSample[1],
        haplotypeParent2 = haplotypesPerSample[2],
        totalHaplotypesInProgeny = sum(haplotypesPerSample[-(parentalIDindexes)]),
        medianHaplotypesPerProgeny = median(haplotypesPerSample[-(parentalIDindexes)]),
        minHaplotypesPerProgeny = min(haplotypesPerSample[-(parentalIDindexes)]),
        maxHaplotypesPerProgeny = max(haplotypesPerSample[-(parentalIDindexes)]),
        totalGenotypesInProgeny = sum(nonMissingGenotypesPerSample[-(parentalIDindexes)]),
        medianGenotypesPerProgeny = median(nonMissingGenotypesPerSample[-(parentalIDindexes)]),
        minGenotypesPerProgeny = min(nonMissingGenotypesPerSample[-(parentalIDindexes)]),
        maxGenotypesPerProgeny = max(nonMissingGenotypesPerSample[-(parentalIDindexes)]),
        totalRecombinations = sum(haplotypesPerSample[-(parentalIDindexes)] - 14),
        medianRecombinationsPerSample = median(haplotypesPerSample[-(parentalIDindexes)] - 14),
        minRecombinationsPerSample = min(haplotypesPerSample[-(parentalIDindexes)] - 14),
        maxRecombinationsPerSample = max(haplotypesPerSample[-(parentalIDindexes)] - 14),
        whichMinRecombinationsPerSample = names(haplotypesPerSample[-(parentalIDindexes)])[which.min(haplotypesPerSample[-(parentalIDindexes)] - 14)],
        whichMaxRecombinationsPerSample = names(haplotypesPerSample[-(parentalIDindexes)])[which.max(haplotypesPerSample[-(parentalIDindexes)] - 14)],
        totalSingleSNPhaplotypes = sum(singleSNPhaplotypesPerSample[-(parentalIDindexes)]),
        medianSingleSNPhaplotypesPerSample = median(singleSNPhaplotypesPerSample[-(parentalIDindexes)]),
        minSingleSNPhaplotypesPerSample = min(singleSNPhaplotypesPerSample[-(parentalIDindexes)]),
        maxSingleSNPhaplotypesPerSample = max(singleSNPhaplotypesPerSample[-(parentalIDindexes)]),
        whichMinSingleSNPhaplotypesPerSample = names(singleSNPhaplotypesPerSample[-(parentalIDindexes)])[which.min(singleSNPhaplotypesPerSample[-(parentalIDindexes)])],
        whichMaxSingleSNPhaplotypesPerSample = names(singleSNPhaplotypesPerSample[-(parentalIDindexes)])[which.max(singleSNPhaplotypesPerSample[-(parentalIDindexes)])],
        totalPutativeCrossovers = sum(putativeCrossoversPerSample[-(parentalIDindexes)]),
        medianPutativeCrossoversPerSample = median(putativeCrossoversPerSample[-(parentalIDindexes)]),
        minPutativeCrossoversPerSample = min(putativeCrossoversPerSample[-(parentalIDindexes)]),
        maxPutativeCrossoversPerSample = max(putativeCrossoversPerSample[-(parentalIDindexes)]),
        whichMinPutativeCrossoversPerSample = names(putativeCrossoversPerSample[-(parentalIDindexes)])[which.min(putativeCrossoversPerSample[-(parentalIDindexes)])],
        whichMaxPutativeCrossoversPerSample = names(putativeCrossoversPerSample[-(parentalIDindexes)])[which.max(putativeCrossoversPerSample[-(parentalIDindexes)])],
        numberOfVariants = dim(vcfFiltered)[1],
        numberOfMendelianErrorVariants = length(which(values(info(vcfFiltered))[["MendelianErrors"]] > 0)),
        numberOfMendelianErrorGenotypes = sum(values(info(vcfFiltered))[["MendelianErrors"]]),
        numberOfSingleSNPhaplotypeVariants = length(which(values(info(vcfFiltered))[["numSingleSNPhaplotypes"]] > 0)),
        numberOfSingleSNPhaplotypesInMultipleSamples = length(which(values(info(vcfFiltered))[["numSingleSNPhaplotypes"]] > 1)),
        meanNumberOfSingleSNPhaplotypesPerSample = sum(values(info(vcfFiltered))[["numSingleSNPhaplotypes"]], na.rm=TRUE) / dim(vcfFiltered)[2],
        numberOfSegregatingSites = length(which(values(info(vcfFiltered))[["SEGREGATING"]])),
        titvRatio = titvRatio,
        titvRatioExcludingAT = titvRatioExcludingAT,
        meanDuplicateDiscordanceRate = meanDuplicateDiscordanceRate,
#        totalRecombinations_Pf3D7_01_v3 = sum(recombinationsPerSample_Pf3D7_01_v3),
#        medianRecombinationsPerSample_Pf3D7_01_v3 = median(recombinationsPerSample_Pf3D7_01_v3),
#        numberOfVariants_Pf3D7_01_v3 = dim(vcf_Pf3D7_01_v3)[1],
#        numberOfMendelianErrors_Pf3D7_01_v3 = length(which(values(info(vcf_Pf3D7_01_v3))[["MendelianErrors"]] > 0)),
#        numberOfSingleSNPhaplotypes_Pf3D7_01_v3 = length(which(values(info(vcf_Pf3D7_01_v3))[["numSingleSNPhaplotypes"]] > 0)),
#        numberOfSegregatingSites_Pf3D7_01_v3 = length(which(values(info(vcf_Pf3D7_01_v3))[["SEGREGATING"]])),
        row.names = basename(theseFiltersPlotFilestem)
#        row.names = paste(c(regionsMaskFilterName, names(additionalInfoFilters)), collapse=".")
      )
      save(returnDF, file=paste(theseFiltersPlotFilestem, "returnDF.rda", sep="."))
    #  save(returnDF, file=paste(paste(c(plotFilestem, regionsMaskFilterName, names(additionalInfoFilters)), collapse="."), "returnDF.rda", sep="."))
      returnDF
    }
  )
  if(shouldReturnVcfOnly) {
    return(resultsList)
  }
  fullReturnDF <- do.call(rbind, resultsList)
  save(
    fullReturnDF,
    file=paste(
      paste(
        c(
          plotFilestem,
          sapply(additionalInfoFilters, function(x) x[["value"]])
        ),
        collapse="."
      ),
      "fullReturnDF.rda", sep="."
    )
  )
  return(fullReturnDF)
}
