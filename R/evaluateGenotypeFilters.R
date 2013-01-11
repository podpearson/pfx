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
  shouldSetNonSegregatingFilt = TRUE,
  plotFilestem                = "evaluateGenotypeFilters",
  GTsToIntMapping             = c("0"=1, "1"=2, "."=0, "./."=0),
  sampleIDcolumn              = "ena_run_accession",
  sampleIDmappingsColumn      = sampleIDcolumn,
  sampleDuplicates            = NULL,
  shouldRecalculateDepthSD    = FALSE,
  shouldCalculateExtraQUAL    = FALSE,
  shouldFilterGenotypes       = TRUE,
  shouldCreateQCFilteringPlots= FALSE,
  shouldCreateRecombPlots     = FALSE,
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
      if(shouldFilterGenotypes) {
        vcfFiltered <- filterGenotypes(
          vcfFiltered,
          genotypeFilters             = genotypeFilters
    #      shouldSetNonSegregatingFilt = shouldSetNonSegregatingFilt,
    #      regionsMask                 = regionsMask,
    #      maxNumFilteredGenotypes     = maxNumFilteredGenotypes
        )
      }
      vcfFiltered <- filterVcf(
        setVcfFilters(
          vcf,
          regionsMask                 = regionsMask,
          additionalInfoFilters       = additionalInfoFilters,
          shouldSetNonSegregatingFilt = shouldSetNonSegregatingFilt
        ),
        keepPASSvariantsOnly = TRUE
#        filtersToRemove = c(regionsMaskFilterName, names(additionalInfoFilters))
      )
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
              plotFilestem=paste(c(plotFilestem, maxNumFilteredGenotypes, genotypeFilters[["LowDP"]][["value"]], genotypeFilters[["HighMAF"]][["value"]], errorVariable), collapse="."),
      #        plotFilestem=paste(c(plotFilestem, regionsMaskFilterName, names(additionalInfoFilters), errorVariable), collapse="."),
              shouldCreateErrorRateBySites=FALSE,
              errorVariable=errorVariable
            )
          }
        )
      }
      mgRecombinations <- recombinationPoints(vcfFiltered, shouldCharacterise=FALSE, GTsToIntMapping=GTsToIntMapping)
      recombinationsPerSample <- rev(
        sapply(
          mgRecombinations[["sampleLevelResults"]],
          function(x) {
            sum(sapply(x, length))
          }
        )
      )
      haplotypesPerSampleAndChromosomeMatrix <- haplotypesPerSample(vcfFiltered)
      haplotypesPerSample <- rowSums(haplotypesPerSampleAndChromosomeMatrix)
      if(shouldCreateRecombPlots) {
        recombinationPlotSeries(
          vcfFiltered,
          plotFilestem                = paste(c(plotFilestem, maxNumFilteredGenotypes, genotypeFilters[["LowDP"]][["value"]], genotypeFilters[["HighMAF"]][["value"]]), collapse="."),
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
        totalHaplotypesInProgeny = sum(haplotypesPerSample[-(1:2)]),
        medianHaplotypesPerProgeny = median(haplotypesPerSample[-(1:2)]),
        minHaplotypesPerProgeny = min(haplotypesPerSample[-(1:2)]),
        maxHaplotypesPerProgeny = max(haplotypesPerSample[-(1:2)]),
        totalRecombinations = sum(recombinationsPerSample),
        medianRecombinationsPerSample = median(recombinationsPerSample),
        minRecombinationsPerSample = min(recombinationsPerSample),
        maxRecombinationsPerSample = max(recombinationsPerSample),
        whichMinRecombinationsPerSample = which.min(recombinationsPerSample),
        whichMaxRecombinationsPerSample = which.max(recombinationsPerSample),
        numberOfVariants = dim(vcfFiltered)[1],
        numberOfMendelianErrors = length(which(values(info(vcfFiltered))[["MendelianErrors"]] > 0)),
        numberOfSingleSNPhaplotypes = length(which(values(info(vcfFiltered))[["numSingleSNPhaplotypes"]] > 0)),
        numberOfSegregatingSites = length(which(values(info(vcfFiltered))[["SEGREGATING"]])),
        titvRatio = titvRatio,
        titvRatioExcludingAT = titvRatioExcludingAT,
#        totalRecombinations_Pf3D7_01_v3 = sum(recombinationsPerSample_Pf3D7_01_v3),
#        medianRecombinationsPerSample_Pf3D7_01_v3 = median(recombinationsPerSample_Pf3D7_01_v3),
#        numberOfVariants_Pf3D7_01_v3 = dim(vcf_Pf3D7_01_v3)[1],
#        numberOfMendelianErrors_Pf3D7_01_v3 = length(which(values(info(vcf_Pf3D7_01_v3))[["MendelianErrors"]] > 0)),
#        numberOfSingleSNPhaplotypes_Pf3D7_01_v3 = length(which(values(info(vcf_Pf3D7_01_v3))[["numSingleSNPhaplotypes"]] > 0)),
#        numberOfSegregatingSites_Pf3D7_01_v3 = length(which(values(info(vcf_Pf3D7_01_v3))[["SEGREGATING"]])),
        row.names = paste(c(regionsMaskFilterName, names(additionalInfoFilters)), collapse=".")
      )
      save(returnDF, file=paste(paste(c(plotFilestem, genotypeFilters[["LowGQ"]][["value"]], genotypeFilters[["LowDP"]][["value"]], genotypeFilters[["HighMAF"]][["value"]]), collapse="."), "returnDF.rda", sep="."))
    #  save(returnDF, file=paste(paste(c(plotFilestem, regionsMaskFilterName, names(additionalInfoFilters)), collapse="."), "returnDF.rda", sep="."))
      returnDF
    }
  )
  fullReturnDF <- do.call(rbind, resultsList)
  return(fullReturnDF)
}
