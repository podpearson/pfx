# pipeline2.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################

# install_github("pfx")
# library("pfx")
# run20121221_3d7_hb3_snps <- pipeline2(parentalStrains=c("ERR019061", "ERR019054"), overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
# run20121221_3d7_hb3_indels <- pipeline2(variantType="indels", parentalStrains=c("ERR019061", "ERR019054"))
# run20121221_7g8_gb4_snps <- pipeline2("7g8_gb4", parentalStrains=c("ERR027099", "ERR027100"))
# run20121221_7g8_gb4_indels <- pipeline2("7g8_gb4", variantType="indels", parentalStrains=c("ERR027099", "ERR027100"))
# run20121221_3d7_hb3_snps <- pipeline2("hb3_dd2", parentalStrains=c("ERR012788", "ERR012840"))
# run20121221_3d7_hb3_indels <- pipeline2("hb3_dd2", variantType="indels", parentalStrains=c("ERR012788", "ERR012840"))

# run20130104_3d7_hb3_snps <- pipeline2(parentalStrains=c("ERR019061", "ERR019054"), genotypesFileFmt="%s.annotated.vcf", overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
# run20130104_3d7_hb3_indels <- pipeline2(variantType="indels", parentalStrains=c("ERR019061", "ERR019054"), genotypesFileFmt="%s.annotated.vcf", overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
# run20130104_7g8_gb4_snps <- pipeline2("7g8_gb4", parentalStrains=c("ERR027099", "ERR027100"), genotypesFileFmt="%s.annotated.vcf", overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
# run20130104_7g8_gb4_indels <- pipeline2("7g8_gb4", variantType="indels", parentalStrains=c("ERR027099", "ERR027100"), genotypesFileFmt="%s.annotated.vcf", overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
# run20130104_hb3_dd2_snps <- pipeline2("hb3_dd2", parentalStrains=c("ERR012788", "ERR012840"), genotypesFileFmt="%s.annotated.vcf", overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
# run20130104_hb3_dd2_indels <- pipeline2("hb3_dd2", variantType="indels", parentalStrains=c("ERR012788", "ERR012840"), genotypesFileFmt="%s.annotated.vcf", overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)

# run20130107_3d7_hb3_snps <- pipeline2(parentalStrains=c("ERR019061", "ERR019054"), genotypesFileFmt="%s.annotated.vcf", overwriteExisting=TRUE, shouldUseSavedVersions=FALSE, genotypesDirectory="data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20120107/gatk")
# run20130107_3d7_hb3_indels <- pipeline2(variantType="indels", parentalStrains=c("ERR019061", "ERR019054"), genotypesFileFmt="%s.annotated.vcf", overwriteExisting=TRUE, shouldUseSavedVersions=FALSE, genotypesDirectory="data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20120107/gatk")
# run20130107_7g8_gb4_snps <- pipeline2("7g8_gb4", parentalStrains=c("ERR027099", "ERR027100"), genotypesFileFmt="%s.annotated.vcf", overwriteExisting=TRUE, shouldUseSavedVersions=FALSE, genotypesDirectory="data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20120107/gatk")
# run20130107_7g8_gb4_indels <- pipeline2("7g8_gb4", variantType="indels", parentalStrains=c("ERR027099", "ERR027100"), genotypesFileFmt="%s.annotated.vcf", overwriteExisting=TRUE, shouldUseSavedVersions=FALSE, genotypesDirectory="data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20120107/gatk")
# run20130107_hb3_dd2_snps <- pipeline2("hb3_dd2", parentalStrains=c("ERR012788", "ERR012840"), genotypesFileFmt="%s.annotated.vcf", overwriteExisting=TRUE, shouldUseSavedVersions=FALSE, genotypesDirectory="data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20120107/gatk")
# run20130107_hb3_dd2_indels <- pipeline2("hb3_dd2", variantType="indels", parentalStrains=c("ERR012788", "ERR012840"), genotypesFileFmt="%s.annotated.vcf", overwriteExisting=TRUE, shouldUseSavedVersions=FALSE, genotypesDirectory="data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20120107/gatk")

pipeline2 <- function(
  cross                       = "3d7_hb3",
  variantType                 = "snps",
  genotypesDirectory          = "data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis/gatk",
  chromosomes                 = sprintf("Pf3D7_%02d_v3", 1:14),
  outputDirectory             = genotypesDirectory,
#  outputDirectory             = sub("genotypes", "genotypes_analysis", genotypesDirectory),
  genotypesFileFmt            = "%s.raw.vcf",
  gffFilename                 = "/data/malariagen2/plasmodium/pf-crosses/data/genome/sanger/version3/September_2012/Pf3D7_v3.gatk.gff",
  gffGRL                      = readGffAsGRangesList(gffFilename, chromsomeNames=chromosomes),
  parentalStrains             = NULL,
  sampleIDcolumn              = "ena_run_accession",
  sampleIDmappingsColumn      = sampleIDcolumn,
  filtersToRemove             = NULL,
  samplesToRemove             = NULL,
#  samplesToRemove             = c("ERR045643", "ERR045644", "ERR045645", "ERR045646", "ERR045647", "ERR045625"),
  plotFilestem                = NULL,
  heterozygosityThreshold     = 0.3,
  missingnessThreshold        = 0.3,
  discordanceThresholdInitial = 1000,
  discordanceThresholdMg      = 10000,
  discordanceThresholdSeg     = 100,
  discordanceThresholdJiang   = 100,
  discordanceThresholdRawVsJia= 100,
  discordanceThresholdFltVsJia= 100,
  discordanceThreshold5ChVsJia= 50,
  shouldCompareWithJiang      = grepl("7g8xGb4", cross),
  overwriteExisting           = NULL,
  shouldUseSavedVersions      = TRUE,
  vcfListRda                  = file.path(outputDirectory, cross, variantType, paste(cross, "vcfList.rda", sep=".")),
  vcfVariantRda               = file.path(outputDirectory, cross, variantType, paste(cross, "vcfVariant.rda", sep=".")),
  vcfVariantAnnotatedRda      = file.path(outputDirectory, cross, variantType, paste(cross, "vcfVariantAnnotated.rda", sep=".")),
  vcfInitialFilteredRda       = file.path(outputDirectory, cross, variantType, paste(cross, "vcfInitialFiltered.rda", sep=".")),
  vcfFinalFilteredRda         = file.path(outputDirectory, cross, variantType, paste(cross, "vcfFinalFiltered.rda", sep=".")),
  vcfUnfilteredFinalSamplesRda = file.path(outputDirectory, cross, variantType, paste(cross, "vcfUnfilteredFinalSamples.rda", sep=".")),
  vcfCoreFinalSamplesRda      = file.path(outputDirectory, cross, variantType, paste(cross, "vcfCoreFinalSamples.rda", sep=".")),
  vcfUniqueFilteredRda        = file.path(outputDirectory, cross, variantType, paste(cross, "vcfUniqueFiltered.rda", sep="."))
) {
  if(file.exists(vcfListRda) & shouldUseSavedVersions) {
    load(vcfListRda)
  } else {
    if(is.null(overwriteExisting)) {
      vcfList <- variantSitesRdaFiles2(
        cross                       = cross,
        variantType                 = variantType,
        genotypesDirectory          = genotypesDirectory,
        chromosomes                 = chromosomes,
        outputDirectory             = outputDirectory,
        genotypesFileFmt            = genotypesFileFmt,
        filtersToRemove             = filtersToRemove,
        samplesToRemove             = samplesToRemove,
        parentalStrains             = parentalStrains
      )
    } else {
      vcfList <- variantSitesRdaFiles2(
        cross                       = cross,
        variantType                 = variantType,
        genotypesDirectory          = genotypesDirectory,
        chromosomes                 = chromosomes,
        outputDirectory             = outputDirectory,
        genotypesFileFmt            = genotypesFileFmt,
        filtersToRemove             = filtersToRemove,
        samplesToRemove             = samplesToRemove,
        overwriteExisting           = overwriteExisting,
        parentalStrains             = parentalStrains
      )
    }
    save(vcfList, file=vcfListRda)
  }
  if(file.exists(vcfVariantRda) & shouldUseSavedVersions) {
    load(vcfVariantRda)
  } else {
    vcfVariant <- combineVcfListIntoVcf(vcfList)
    save(vcfVariant, file=vcfVariantRda)
  }
  if(file.exists(vcfVariantAnnotatedRda) & shouldUseSavedVersions) {
    load(vcfVariantAnnotatedRda)
  } else {
    vcfVariantAnnotated <- annotateVcf(vcfVariant)
    save(vcfVariantAnnotated, file=vcfVariantAnnotatedRda)
  }
  if(file.exists(vcfInitialFilteredRda) & shouldUseSavedVersions) {
    load(vcfInitialFilteredRda)
  } else {
    vcfInitialFiltered <- setVcfFilters(
      vcfVariantAnnotated,
      regionsMask                 = varRegions_v3(),
      additionalInfoFilters = list(
        "LowQD" = list(column="QD", operator="<=", value=36),
        "HighSB" = list(column="SB", operator=">=", value=-6000)
      )
#      shouldSetMultiallelicFilter = TRUE,
#      shouldSetNonSegregatingFilt = TRUE
    )
    save(vcfInitialFiltered, file=vcfInitialFilteredRda)
  }
  vcfSegregating <- filterVcf(vcfInitialFiltered, keepPASSvariantsOnly=TRUE)
  initialSampleQCresults <- sampleQC(
    vcfSegregating,
    discordanceThreshold=discordanceThresholdInitial,
    plotFilestem=file.path(outputDirectory, cross, variantType, paste(cross, "initital", sep=".")),
    gffGRL=gffGRL,
    sampleIDcolumn=sampleIDcolumn,
    sampleIDmappingsColumn=sampleIDmappingsColumn
  )
  initialSNPnumbersMatrix <- recombinationPlotSeries(
    vcfInitialFiltered,
    plotFilestem=file.path(outputDirectory, cross, variantType, paste(cross, "allSamples", sep=".")),
    filters=c("InVarRegion", "LowQD", "HighSB"),
    sampleIDcolumn=sampleIDcolumn,
    sampleIDmappingsColumn=sampleIDmappingsColumn,
    sampleDuplicates=initialSampleQCresults[["sampleDuplicates"]]
  )
  finalSamples <- setdiff(dimnames(vcfInitialFiltered)[[2]], initialSampleQCresults[["qcFailedSamples"]])
#  vcfUnfilteredFinalSamples <- annotateVcf(vcfVariant[, finalSamples])
#  save(vcfUnfilteredFinalSamples, file=vcfUnfilteredFinalSamplesRda)
  vcfCoreFinalSamples <- filterVcf(
    setVcfFilters(
      annotateVcf(vcfVariant[, finalSamples]),
      regionsMask                 = varRegions_v3(),
    ),
    filtersToRemove = "InVarRegion"
  )
  save(vcfCoreFinalSamples, file=vcfCoreFinalSamplesRda)
  coreVcfFinalSamples <- annotateVcf(
    filterVcf(
      setVcfFilters(
        vcfVariantAnnotated[, finalSamples],
        regionsMask                 = varRegions_v3(),
        additionalInfoFilters = list(
          "LowQD" = list(column="QD", operator="<=", value=36),
          "HighSB" = list(column="SB", operator=">=", value=-6000)
        )
  #      shouldSetMultiallelicFilter = TRUE,
  #      shouldSetNonSegregatingFilt = TRUE
      ),
      filtersToRemove = "InVarRegion"
    )
  )
  qcFilteringResults_coreFinalSamples <- qcFilteringPlots(coreVcfFinalSamples, plotFilestem=file.path(outputDirectory, cross, variantType, paste(cross, "coreFinalSamples", sep=".")))
  qcFilteringResults_coreFinalSamplesMaxMAF <- qcFilteringPlots(coreVcfFinalSamples, plotFilestem=file.path(outputDirectory, cross, variantType, paste(cross, "coreFinalSamplesMaxMAF", sep=".")), errorVariable="maxMAF", errorThreshold=0.1)
  finalSNPnumbersMatrix <- recombinationPlotSeries(
    coreVcfFinalSamples,
    plotFilestem=file.path(outputDirectory, cross, variantType, paste(cross, "coreFinalSamples", sep=".")),
    filters=c("InVarRegion", "LowQD", "HighSB"),
    sampleIDcolumn=sampleIDcolumn,
    sampleIDmappingsColumn=sampleIDmappingsColumn,
    sampleDuplicates=initialSampleQCresults[["sampleDuplicates"]]
  )
  vcfSegregating <- filterVcf(coreVcfFinalSamples, keepPASSvariantsOnly=TRUE)
    
  gc()
  if(shouldCompareWithJiang) {
    jiangVcf <- loadJiangGenotypesAsVcf() # make this return a VCF?
    jiangSampleQCresults <- sampleQC(jiangVcf, discordanceThresholdJiang, plotFilestem="JiangEtAl", GTsToIntMapping = c("7"=1, "G"=2, "."=0), shouldCalcMissingnessAndHet=FALSE, shouldRenameSamples=FALSE) # should output heatmap of discordances
#    jiangUniqueSamples <- uniqueSamples(jiangVcf, discordanceThresholdJiang, plotFilestem="JiangEtAl", GTsToIntMapping = c("7"=1, "G"=2, "."=0)) # should output heatmap of discordances
    if(cross=="v3UG_7g8xGb4") {
      seqlevels(jiangVcf) <- sprintf("Pf3D7_%02d_v3", as.integer(sub("MAL", "", seqlevels(jiangVcf))))
    }
    genotypeConcordanceRaw <- compareCalls(vcfVariantAnnotated, jiangVcf, plotFilestem=file.path(outputDirectory, cross, variantType, paste(cross, "comparison", "raw", sep=".")), discordanceThreshold=discordanceThresholdRawVsJia)
#    genotypeConcordanceRaw <- compareCalls(vcfFiltered, jiangVcf, plotFilestem=paste(cross, "comparison", "raw", sep="."), discordanceThreshold=discordanceThresholdRawVsJia)
    genotypeConcordance <- compareCalls(coreVcfFinal, jiangVcf, plotFilestem=file.path(outputDirectory, cross, variantType, paste(cross, "comparison", "filtered", sep=".")), discordanceThreshold=discordanceThresholdFltVsJia) # Should give slide 3, histogram of pair-wise numbers of discordant, heatmap of sample discordances and heatmap for discordances for presumed identical, recombinationPlot of both together
#    genotypeConcordance <- compareCalls(vcfSegregating, jiangVcf, plotFilestem=paste(cross, "comparison", "filtered", sep="."), discordanceThreshold=discordanceThresholdFltVsJia) # Should give slide 3, histogram of pair-wise numbers of discordant, heatmap of sample discordances and heatmap for discordances for presumed identical, recombinationPlot of both together
#    genotypeConcordancePf3D7_02_v3 <- compareCalls(vcfSegregating[seqnames(vcfSegregating)=="Pf3D7_02_v3"], jiangVcf[seqnames(jiangVcf)=="Pf3D7_02_v3"], plotFilestem=paste(cross, "comparison", "Pf3D7_02_v3", sep="."), discordanceThreshold=discordanceThreshold5ChVsJia) # Should give slide 3, histogram of pair-wise numbers of discordant, heatmap of sample discordances and heatmap for discordances for presumed identical, recombinationPlot of both together
#    genotypeConcordance6chromosomes <- compareCalls(vcfSegregating[seqnames(vcfSegregating) %in% c("Pf3D7_02_v3", "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3")], jiangVcf[seqnames(jiangVcf) %in% c("Pf3D7_02_v3", "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3")], plotFilestem=paste(cross, "comparison", "6chromosomes", sep="."), discordanceThreshold=discordanceThreshold5ChVsJia) # Should give slide 3, histogram of pair-wise numbers of discordant, heatmap of sample discordances and heatmap for discordances for presumed identical, recombinationPlot of both together
#    genotypeConcordance5chromosomes <- compareCalls(vcfSegregating[seqnames(vcfSegregating) %in% c("Pf3D7_02_v3", "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_10_v3", "Pf3D7_11_v3")], jiangVcf[seqnames(jiangVcf) %in% c("Pf3D7_02_v3", "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_10_v3", "Pf3D7_11_v3")], plotFilestem=paste(cross, "comparison", "5chromosomes", sep="."), discordanceThreshold=discordanceThreshold5ChVsJia) # Should give slide 3, histogram of pair-wise numbers of discordant, heatmap of sample discordances and heatmap for discordances for presumed identical, recombinationPlot of both together
    gc()
    if(!file.exists(file.path(outputDirectory, cross, variantType, paste(cross, "mgRecombinations.rda", sep="."))) | !shouldUseSavedVersions) {
      mgRecombinations <- recombinationPoints(vcfSegregating[filt(vcfSegregating)=="PASS", qcPlusUniqueSamples], gffGRL) # extend crossoversAnalysis to include classification as exonic, intronic, etc
      save(mgRecombinations, file=file.path(outputDirectory, cross, variantType, paste(cross, "mgRecombinations.rda", sep=".")))
    } else {
      load(file.path(outputDirectory, cross, variantType, paste(cross, "mgRecombinations.rda", sep=".")))
    }
    if(!file.exists("~/jiangRecombinations.rda") | !shouldUseSavedVersions) {
      jiangRecombinations <- recombinationPoints(jiangVcf, gffGRL, GTsToIntMapping = c("7"=1, "G"=2, "."=0))
      save(jiangRecombinations, file="jiangRecombinations.rda")
    } else {
      load("~/jiangRecombinations.rda")
    }
    medianBreakpointAccuracies <- compareRecombinations(mgRecombinations, jiangRecombinations) # to include slide 9 plot, slide 12 plot, median accuracies, Venn
  } else {
    if(!file.exists(file.path(outputDirectory, cross, variantType, paste(cross, "mgRecombinations.rda", sep="."))) | !shouldUseSavedVersions) {
      mgRecombinations <- recombinationPoints(vcfSegregating[filt(vcfSegregating)=="PASS"], gffGRL) # extend crossoversAnalysis to include classification as exonic, intronic, etc
      save(mgRecombinations, file=file.path(outputDirectory, cross, variantType, paste(cross, "mgRecombinations.rda", sep=".")))
    } else {
      load(file.path(outputDirectory, cross, variantType, paste(cross, "mgRecombinations.rda", sep=".")))
    }
  }
  recombinationRates <- analyseRecombinations(mgRecombinations, plotFilestem=file.path(outputDirectory, cross, variantType, cross)) # to include slide 13 plot, breakdown of CO and GC by progeny, chromosome and by cross, CO and GC rates
  returnList <- list(
    vcfVariantAnnotated                  = vcfVariantAnnotated,
    vcfSegregating                       = vcfSegregating,
    mgRecombinations                     = mgRecombinations,
    recombinationRates                   = recombinationRates,
    qcFilteringResults_coreFinalSamples  = qcFilteringResults_coreFinalSamples,
    qcFilteringResults_coreFinalSamplesMaxMAF = qcFilteringResults_coreFinalSamplesMaxMAF,
#    qcFilteringResults                   = qcFilteringResults,
    initialSampleQCresults               = initialSampleQCresults,
    initialSNPnumbersMatrix              = initialSNPnumbersMatrix,
#    qcFilteringResultsFinalPostFiltering = qcFilteringResultsFinalPostFiltering,
#    finalSampleQCresults                 = finalSampleQCresults,
    finalSNPnumbersMatrix                = finalSNPnumbersMatrix
#    uniqueSNPnumbersMatrix               = uniqueSNPnumbersMatrix,
#    finalUniqueSampleQCresults           = finalUniqueSampleQCresults
  )
  if(shouldCompareWithJiang) {
    returnList <- c(
      returnList,
      list(
        jiangSampleQCresults                 = jiangSampleQCresults,
        jiangRecombinations                  = jiangRecombinations,
        genotypeConcordance                  = genotypeConcordance,
#        genotypeConcordancePf3D7_02_v3       = genotypeConcordancePf3D7_02_v3,
        medianBreakpointAccuracies           = medianBreakpointAccuracies
      )
    )
  }
  save(returnList, file=file.path(outputDirectory, cross, variantType, paste(cross, "returnList.rda", sep=".")))
  writeVcf(vcfSegregating, filename=file.path(outputDirectory, cross, variantType, paste(cross, "filtered.vcf", sep=".")), index=TRUE)
  return(returnList)
}

