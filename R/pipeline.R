# pipeline.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################

#temp_7g8xGb4 <- pipeline(file.path("/data/galton/users/rpearson", "crossesTesting", "release", "7g8xGb4-qcPlusSamples-0.1.vcf.gz"))
#temp_3d7xHb3 <- pipeline(file.path("/data/galton/users/rpearson", "crossesTesting", "release", "3d7xHb3-qcPlusSamples-0.1.vcf.gz"), discordanceThresholdMg=5000)
#temp_Hb3xDd2 <- pipeline(file.path("/data/galton/users/rpearson", "crossesTesting", "release", "Hb3xDd2-qcPlusSamples-0.1.vcf.gz"))
#pipeline_7g8xGb4 <- pipeline("7g8xGb4")
#pipeline_3d7xHb3 <- pipeline("3d7xHb3", discordanceThresholdMg=5000)
#pipeline_Hb3xDd2 <- pipeline("Hb3xDd2")
#v3UG_7g8xGb4 <- pipeline("v3UG_7g8xGb4", "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_default/snp_genotypes_analysis/PFproj3-unigen-snponly-BQ20.vcf.gz", chromosomes=sprintf("Pf3D7_%02d_v3", 1:14), filtersToRemove="LowQual", samplesToRemove = c("ERR045643", "ERR045644", "ERR045645", "ERR045646", "ERR045647", "ERR045625"), overwriteExisting=TRUE, discordanceThresholdMg=28000, plotFilestem="v3UG_7g8xGb4")
#v3UG_7g8xGb4 <- pipeline("v3UG_7g8xGb4", "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_default/snp_genotypes_analysis/PFproj3-unigen-snponly-BQ20.vcf.gz", chromosomes=sprintf("Pf3D7_%02d_v3", 1:14), filtersToRemove="LowQual", samplesToRemove = c("ERR045643", "ERR045644", "ERR045645", "ERR045646", "ERR045647", "ERR045625"), overwriteExisting=TRUE, discordanceThresholdMg=28000, plotFilestem="v3UG_7g8xGb4", discordanceThresholdRawVsJia=300)
#v3UG_7g8xGb4 <- pipeline("v3UG_7g8xGb4", "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_default/snp_genotypes_analysis/PFproj3-unigen-snponly-BQ20.vcf.gz", chromosomes=sprintf("Pf3D7_%02d_v3", 1:14), filtersToRemove="LowQual", samplesToRemove = c("ERR045643", "ERR045644", "ERR045645", "ERR045646", "ERR045647", "ERR045625", "ERR029145"), overwriteExisting=TRUE, discordanceThresholdMg=28000, plotFilestem="v3UG_7g8xGb4", discordanceThresholdRawVsJia=300, shouldUseSavedVersions=FALSE)
#v3UG_7g8xGb4 <- pipeline("v3UG_7g8xGb4", "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_default/snp_genotypes_analysis/PFproj3-unigen-snponly-BQ20.vcf.gz", chromosomes=sprintf("Pf3D7_%02d_v3", 1:14), filtersToRemove="LowQual", samplesToRemove = c("ERR045643", "ERR045644", "ERR045645", "ERR045646", "ERR045647", "ERR045625", "ERR029145"), overwriteExisting=TRUE, discordanceThresholdMg=28000, plotFilestem="v3UG_7g8xGb4", discordanceThresholdRawVsJia=300)

#v3UG_7g8xGb4 <- pipeline("v3UG_7g8xGb4", "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_default/snp_genotypes_analysis/PFproj3-unigen-snponly-BQ20.vcf.gz", chromosomes=sprintf("Pf3D7_%02d_v3", 1:14), filtersToRemove="LowQual", samplesToRemove = c("ERR045643", "ERR045644", "ERR045645", "ERR045646", "ERR045647", "ERR045625"), discordanceThresholdMg=28000, plotFilestem="v3UG_7g8xGb4", discordanceThresholdRawVsJia=300, overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
#v3UG_3d7xHb3 <- pipeline("v3UG_3d7xHb3", "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_default/snp_genotypes_analysis/PFproj2-unigen-snponly-BQ20.vcf.gz", chromosomes=sprintf("Pf3D7_%02d_v3", 1:14), filtersToRemove="LowQual", discordanceThresholdMg=15000, plotFilestem="v3UG_3d7xHb3", parentalStrains=c("ERR019061", "ERR019054"), overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)
#v3UG_Hb3xDd2 <- pipeline("v3UG_Hb3xDd2", "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_default/snp_genotypes_analysis/PFproj1-unigen-snponly-BQ20.vcf.gz", chromosomes=sprintf("Pf3D7_%02d_v3", 1:14), filtersToRemove="LowQual", discordanceThresholdMg=28000, plotFilestem="v3UG_Hb3xDd2", overwriteExisting=TRUE, shouldUseSavedVersions=FALSE)

#v3UG_3d7xHb3 <- pipeline("v3UG_3d7xHb3", "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_default/snp_genotypes_analysis/PFproj2-unigen-snponly-BQ20.vcf.gz", chromosomes=sprintf("Pf3D7_%02d_v3", 1:14), filtersToRemove="LowQual", discordanceThresholdMg=15000, plotFilestem="v3UG_3d7xHb3", parentalStrains=c("ERR019061", "ERR019054"), overwriteExisting=FALSE, shouldUseSavedVersions=TRUE)
#v3UG_Hb3xDd2 <- pipeline("v3UG_Hb3xDd2", "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_default/snp_genotypes_analysis/PFproj1-unigen-snponly-BQ20.vcf.gz", chromosomes=sprintf("Pf3D7_%02d_v3", 1:14), filtersToRemove="LowQual", discordanceThresholdMg=28000, plotFilestem="v3UG_Hb3xDd2", overwriteExisting=FALSE, shouldUseSavedVersions=TRUE)
#v3UG_7g8xGb4 <- pipeline("v3UG_7g8xGb4", "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_default/snp_genotypes_analysis/PFproj3-unigen-snponly-BQ20.vcf.gz", chromosomes=sprintf("Pf3D7_%02d_v3", 1:14), filtersToRemove="LowQual", samplesToRemove = c("ERR045643", "ERR045644", "ERR045645", "ERR045646", "ERR045647", "ERR045625"), discordanceThresholdMg=28000, plotFilestem="v3UG_7g8xGb4", discordanceThresholdRawVsJia=300, overwriteExisting=FALSE, shouldUseSavedVersions=TRUE)

pipeline <- function(
  cross                       = "7g8xGb4",
  vcfFilename                 = file.path("/data/galton/users/rpearson", "crossesTesting", "release", paste(cross, "-qcPlusSamples-0.1.vcf.gz", sep="")),
#  chromosomes                 = sprintf("MAL%d", 1:14),
  chromosomes                 = sprintf("Pf3D7_%02d_v3", 1:14),
#  gffFilename                 = "/data/galton/mirror/nfs/team112/annotation/plasmodium/falciparum/Pfalciparum_PlasmoDB-7.2.gff",
  gffFilename                 = "/data/malariagen2/plasmodium/pf-crosses/data/genome/sanger/version3/September_2012/Pf3D7_v3.gatk.gff",
  gffGRL                      = readGffAsGRangesList(gffFilename, chromsomeNames=chromosomes),
  parentalStrains             = NULL,
  sampleIDcolumn              = "ena_run_accession",
  sampleIDmappingsColumn      = sampleIDcolumn,
#  chromosomes                 = sprintf("Pf3D7_%02d_v3", 1:14),
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
  vcfListRda                  = paste(cross, "vcfList.rda", sep="."),
  vcfVariantRda               = paste(cross, "vcfVariant.rda", sep="."),
  vcfVariantAnnotatedRda      = paste(cross, "vcfVariantAnnotated.rda", sep="."),
  vcfInitialFilteredRda       = paste(cross, "vcfInitialFiltered.rda", sep="."),
  vcfFinalFilteredRda         = paste(cross, "vcfFinalFiltered.rda", sep="."),
  vcfUniqueFilteredRda        = paste(cross, "vcfUniqueFiltered.rda", sep=".")
) {
  if(file.exists(vcfListRda) & shouldUseSavedVersions) {
    load(vcfListRda)
  } else {
    if(is.null(overwriteExisting)) {
      vcfList <- variantSitesRdaFiles(vcfFilename, chromosomes=chromosomes, filtersToRemove=filtersToRemove, samplesToRemove=samplesToRemove, parentalStrains=parentalStrains)
    } else {
      vcfList <- variantSitesRdaFiles(vcfFilename, chromosomes=chromosomes, filtersToRemove=filtersToRemove, samplesToRemove=samplesToRemove, overwriteExisting=overwriteExisting, parentalStrains=parentalStrains)
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
  qcFilteringResults_raw <- qcFilteringPlots(vcfVariantAnnotated, plotFilestem=paste(cross, "allSamples", sep="."))
  if(file.exists(vcfInitialFilteredRda) & shouldUseSavedVersions) {
    load(vcfInitialFilteredRda)
  } else {
    vcfInitialFiltered <- setVcfFilters(
      vcfVariantAnnotated,
      regionsMask                 = varRegions_v3(),
      additionalInfoFilters = list(
#        "LowQD" = list(column="QD", operator="<=", value=36),
        "HighMaxMAF" = list(column="maxMAF", operator=">=", value=0.1, filterOutNAs=TRUE),
#        "HighMissingness" = list(column="missingness", operator=">=", value=1, filterOutNAs=TRUE),
        "LowDepth" = list(column="missingness2", operator=">=", value=1, filterOutNAs=TRUE)
      ),
      shouldSetMultiallelicFilter = TRUE,
      shouldSetNonSegregatingFilt = TRUE
    )
    save(vcfInitialFiltered, file=vcfInitialFilteredRda)
  }
  vcfSegregating <- filterVcf(vcfInitialFiltered, keepPASSvariantsOnly=TRUE)
  initialSampleQCresults <- sampleQC(
    vcfSegregating,
    discordanceThreshold=discordanceThresholdInitial,
    plotFilestem=paste(cross, "initital", sep="."),
    gffGRL=gffGRL,
    sampleIDcolumn=sampleIDcolumn,
    sampleIDmappingsColumn=sampleIDmappingsColumn
  )
  initialSNPnumbersMatrix <- recombinationPlotSeries(
    vcfInitialFiltered,
    plotFilestem=paste(cross, "allSamples", sep="."),
#    filters=c("InVarRegion", "LowQD", "LowDepth")
    filters=c("InVarRegion", "HighMaxMAF", "LowDepth", "MultiAllelic", "NonSegregating"),
#    filters=c("InVarRegion", "HighMaxMAF", "LowDepth", "HighMissingness")
    sampleIDcolumn=sampleIDcolumn,
    sampleIDmappingsColumn=sampleIDmappingsColumn,
    sampleDuplicates=initialSampleQCresults[["sampleDuplicates"]]
  )
  coreVcf <- filterVcf(vcfInitialFiltered, filtersToRemove = "InVarRegion")
  qcFilteringResults_core <- qcFilteringPlots(coreVcf, plotFilestem=paste(cross, "core", sep="."))
#  naiveFilteringVcf <- filterVcf(vcfInitialFiltered, filtersToRemove = "InVarRegion", "HighMaxMAF")
#  qcFilteringResults_core <- qcFilteringPlots(coreVcf, plotFilestem=paste(cross, "core", sep="."))
  
#  if(file.exists(vcfFinalFilteredRda) & shouldUseSavedVersions) {
#    load(vcfFinalFilteredRda)
#  } else {
#    if(length(initialSampleQCresults[["qcFailedSamples"]]) > 0) {
#      finalSamples <- setdiff(dimnames(vcfInitialFiltered)[[2]], initialSampleQCresults[["qcFailedSamples"]])
#      vcfVariantOnFinalSamples <- annotateVcf(filterVcf(vcfInitialFiltered[, finalSamples]))
#      vcfFinalSamples <- annotateVcf(vcfVariantOnFinalSamples[, finalSamples])
#    } else {
#      vcfFinalSamples <- vcfInitialFiltered
#    }
#    qcFilteringResultsFinalPreFiltering <- qcFilteringPlots(vcfFinalSamples, plotFilestem=paste(cross, "finalPreFiltering", sep="."))
#    filt(vcfFinalSamples) <- "PASS"
#    vcfFinalFiltered <- setVcfFilters(
#      vcfFinalSamples,
##      additionalInfoFilters = NULL,
#      additionalInfoFilters = list(
##        "LowQD" = list(column="QD", operator="<=", value=1),
#        "HighSB" = list(column="SB", operator=">=", value=-4200),
#        "HighMeanMAF" = list(column="meanMAF", operator=">=", value=0.1, filterOutNAs=TRUE),
#        "HighMaxMAF" = list(column="maxMAF", operator=">=", value=0.1, filterOutNAs=TRUE),
#        "HighMaxParentMAF" = list(column="maxParentMAF", operator=">=", value=0.1, filterOutNAs=TRUE),
#        "HighMissingness" = list(column="missingness2", operator=">=", value=1, filterOutNAs=TRUE),
#        "HighMQ0" = list(column="MQ0", operator=">=", value=10000, filterOutNAs=TRUE)
#      ),
#      regionsMask                 = varRegions_v3(),
#      shouldSetMultiallelicFilter = TRUE,
#      shouldSetNonSegregatingFilt = TRUE,
#      shouldSetMaxNoCallsFilter   = TRUE
#    )
#    save(vcfFinalFiltered, file=vcfFinalFilteredRda)
#  }
#  nonVarVcf <- filterVcf(vcfFinalFiltered, filtersToRemove = "InVarRegion")
#  qcFilteringResultsHypothesisDriven <- qcFilteringPlots(
#    nonVarVcf,
##    filterVcf(vcfFinalFiltered, filtersToRemove = "InVarRegion"),
#    variablesToPlot             = c(
#      "BaseQRankSum"   = "highIsGood",
##      "DS"             = "lowIsGood",
##      "Dels"           = "lowIsGood",
#      "FS"             = "lowIsGood",
#      "HaplotypeScore" = "lowIsGood",
#      "MQ"             = "highIsGood",
#      "MQ0"            = "lowIsGood",
#      "MQRankSum"      = "highIsGood",
#      "QD"             = "highIsGood",
##      "RPA"            = "lowIsGood",
#      "ReadPosRankSum" = "highIsGood",
#      "SB"             = "lowIsGood",
##      "homopolymer5Proximity" = "highIsGood",
#      "homopolymer10Proximity" = "highIsGood",
##      "homopolymer15Proximity" = "highIsGood",
#      "UQ"                     = "lowIsGood",
#      "GC500"                  = "highIsGood"
#  #    "meanMAF"        = "lowIsGood",
#  #    "maxMAF"         = "lowIsGood",
#  #    "maxParentMAF"   = "lowIsGood",
#  #    "missingness"    = "lowIsGood",
#  #    "missingness2"   = "lowIsGood",
#  #    "heterozgosity"  = "lowIsGood"
#    ),
#    errorVariable="maxMAF",
#    errorThreshold=0.1,
#    ylim                        = c(-0.2,0),
#    plotFilestem=paste(cross, "hypothesesDriven", sep=".")
#  )
#  nonVarSampleQCResults <- sampleQC(nonVarVcf, discordanceThreshold=discordanceThresholdInitial, plotFilestem=paste(cross, "nonVar", sep="."), gffGRL=gffGRL)
#  nonVarSNPnumbersMatrix <- recombinationPlotSeries(
#    nonVarVcf,
#    plotFilestem=paste(cross, "nonVar", sep="."),
#    filters=c("InVarRegion", "HighMaxMAF"),
#    sampleIDcolumn=sampleIDcolumn,
#    sampleIDmappingsColumn=sampleIDmappingsColumn,
#    sampleDuplicates=nonVarSampleQCResults[["sampleDuplicates"]]
#  )
#  nonVarVcfUQ40 <- setVcfFilters(nonVarVcf, additionalInfoFilters=list("HighUQ"= list(column="UQ", operator=">=", value=40)))
#  nonVarSNPnumbersMatrix <- recombinationPlotSeries(
#    nonVarVcfUQ40,
#    plotFilestem=paste(cross, "nonVarVcfUQ40", sep="."),
#    filters=c("InVarRegion", "HighUQ"),
#    sampleIDcolumn=sampleIDcolumn,
#    sampleIDmappingsColumn=sampleIDmappingsColumn,
#    sampleDuplicates=nonVarSampleQCResults[["sampleDuplicates"]]
#  )
#  qcFilteringResultsHypothesisDriven <- qcFilteringPlots(
#    filterVcf(nonVarVcfUQ40, filtersToRemove = c("InVarRegion", "HighUQ")),
##    filterVcf(vcfFinalFiltered, filtersToRemove = "InVarRegion"),
#    variablesToPlot             = c(
#      "BaseQRankSum"   = "highIsGood",
#      "FS"             = "lowIsGood",
#      "HaplotypeScore" = "lowIsGood",
#      "MQ"             = "highIsGood",
#      "MQ0"            = "lowIsGood",
#      "MQRankSum"      = "highIsGood",
#      "QD"             = "highIsGood",
#      "ReadPosRankSum" = "highIsGood",
#      "SB"             = "lowIsGood",
##      "homopolymer5Proximity" = "highIsGood",
#      "homopolymer10Proximity" = "highIsGood",
##      "homopolymer15Proximity" = "highIsGood",
#      "UQ"                     = "lowIsGood",
#      "GC500"                  = "highIsGood"
#    ),
#    errorVariable="maxMAF",
#    errorThreshold=0.1,
#    ylim                        = c(-4,0),
#    plotFilestem=paste(cross, "nonVarVcfUQ40", sep=".")
#  )
#  nonVarVcfUQ22MQ00 <- setVcfFilters(
#    nonVarVcf,
#    additionalInfoFilters=list(
#      "HighUQ"= list(column="UQ", operator=">=", value=22),
#      "MQ00"= list(column="MQ0", operator=">=", value=1)
#    )
#  )
#  qcFiltering_nonVarVcfUQ22MQ00 <- qcFilteringPlots(
#    filterVcf(nonVarVcfUQ22MQ00, filtersToRemove = c("InVarRegion", "HighUQ", "MQ00")),
##    filterVcf(vcfFinalFiltered, filtersToRemove = "InVarRegion"),
#    variablesToPlot             = c(
#      "BaseQRankSum"   = "highIsGood",
#      "FS"             = "lowIsGood",
#      "HaplotypeScore" = "lowIsGood",
#      "MQ"             = "highIsGood",
#      "MQ0"            = "lowIsGood",
#      "MQRankSum"      = "highIsGood",
#      "QD"             = "highIsGood",
#      "ReadPosRankSum" = "highIsGood",
#      "SB"             = "lowIsGood",
##      "homopolymer5Proximity" = "highIsGood",
#      "homopolymer10Proximity" = "highIsGood",
##      "homopolymer15Proximity" = "highIsGood",
#      "UQ"                     = "lowIsGood",
#      "GC500"                  = "highIsGood"
#    ),
#    errorVariable="maxMAF",
#    errorThreshold=0.1,
#    ylim                        = c(-4,0),
#    plotFilestem=paste(cross, "nonVarVcfUQ22MQ00", sep=".")
#  )
#  nonVarVcfHP10FS20SBUQ22MQ00QD <- setVcfFilters(
#    nonVarVcf,
#    additionalInfoFilters=list(
#      "homopolymer10Proximity"= list(column="homopolymer10Proximity", operator="<=", value=50),
#      "HighFS"= list(column="FS", operator=">=", value=20),
#      "HighSB"= list(column="SB", operator=">=", value=-22000),
#      "HighUQ"= list(column="UQ", operator=">=", value=22),
#      "MQ00"= list(column="MQ0", operator=">=", value=1),
#      "LowQD"= list(column="QD", operator="<=", value=35)
#    )
#  )
#  qcFiltering_nonVarVcfHP10FS20SBUQ22MQ00QD <- qcFilteringPlots(
#    filterVcf(nonVarVcfHP10FS20SBUQ22MQ00QD, filtersToRemove = c("InVarRegion", "homopolymer10Proximity", "HighFS", "HighSB", "HighUQ", "MQ00", "LowQD")),
##    filterVcf(vcfFinalFiltered, filtersToRemove = "InVarRegion"),
#    variablesToPlot             = c(
#      "BaseQRankSum"   = "highIsGood",
#      "FS"             = "lowIsGood",
#      "HaplotypeScore" = "lowIsGood",
#      "MQ"             = "highIsGood",
#      "MQ0"            = "lowIsGood",
#      "MQRankSum"      = "highIsGood",
#      "QD"             = "highIsGood",
##      "ReadPosRankSum" = "highIsGood",
#      "SB"             = "lowIsGood",
#      "homopolymer5Proximity" = "highIsGood",
#      "homopolymer10Proximity" = "highIsGood",
##      "homopolymer15Proximity" = "highIsGood",
#      "UQ"                     = "lowIsGood",
#      "GC500"                  = "highIsGood"
#    ),
#    errorVariable="maxMAF",
#    errorThreshold=0.1,
#    ylim                        = c(-4,0),
#    plotFilestem=paste(cross, "nonVarVcfHP10FS20UQ22MQ00", sep=".")
#  )
#  nonVarVcfHP10FS20SBUQ22MQ00QD_SNPnumbersMatrix <- recombinationPlotSeries(
#    nonVarVcfHP10FS20SBUQ22MQ00QD,
#    plotFilestem=paste(cross, "nonVarVcfHP10UQ22MQ00", sep="."),
#    filters=c("InVarRegion", "homopolymer10Proximity", "HighFS", "HighSB", "HighUQ", "MQ00", "LowQD"),
#    sampleIDcolumn=sampleIDcolumn,
#    sampleIDmappingsColumn=sampleIDmappingsColumn,
#    sampleDuplicates=nonVarSampleQCResults[["sampleDuplicates"]]
#  )
#
#
#
#  qcFilteringResultsFinalPostFiltering <- qcFilteringPlots(
#    filterVcf(
#      vcfFinalFiltered,
#      shouldRemoveInvariant=FALSE,
#      filtersToRemove=c("LowQD", "HighSB", "HighMeanMAF", "HighMissingness", "HighMQ0", "InVarRegion", "ExcessiveNoCalls")
#    ),
#    plotFilestem=paste(cross, "finalPostFiltering", sep=".")
#  )
##  info(vcfFinalFiltered[filt(vcfFinalFiltered) %in% c("NonSegregating")])
##  stem(values(info(vcfFinalFiltered[filt(vcfFinalFiltered) %in% c("NonSegregating", "PASS")]))[["MQ0"]])
##  stem(values(info(vcfFinalFiltered[filt(vcfFinalFiltered) %in% c("NonSegregating", "PASS")]))[["SB"]][values(info(vcfFinalFiltered[filt(vcfFinalFiltered) %in% c("NonSegregating", "PASS")]))[["SB"]]>-Inf])
##  stem(values(info(vcfFinalFiltered[filt(vcfFinalFiltered) %in% c("NonSegregating", "PASS")]))[["SB"]][values(info(vcfFinalFiltered[filt(vcfFinalFiltered) %in% c("NonSegregating", "PASS")]))[["SB"]]>-5000])
##  info(vcfFinalFiltered[filt(vcfFinalFiltered) %in% c("NonSegregating", "PASS")])
#  finalSampleQCresults <- sampleQC(vcfFinalFiltered, discordanceThreshold=discordanceThresholdInitial, plotFilestem=paste(cross, "final", sep="."), gffGRL=gffGRL)
##  finalSNPnumbersMatrix <- recombinationPlotSeries(vcfFinalFiltered, plotFilestem=paste(cross, "final", sep="."), filters=c("LowQD", "InVarRegion", "NonSegregating", "ExcessiveNoCalls"))
##  finalSNPnumbersMatrix <- recombinationPlotSeries(vcfFinalFiltered, plotFilestem=paste(cross, "final", sep="."), filters=c("LowQD", "HighSB", "InVarRegion", "NonSegregating", "ExcessiveNoCalls"))
##  finalSNPnumbersMatrix <- recombinationPlotSeries(vcfFinalFiltered, plotFilestem=paste(cross, "final", sep="."), filters=c("LowQD", "HighSB", "HighMeanMAF", "HighMissingness", "InVarRegion", "NonSegregating", "ExcessiveNoCalls"))
##  finalSNPnumbersMatrix <- recombinationPlotSeries(vcfFinalFiltered, plotFilestem=paste(cross, "final", sep="."), filters=c("HighMissingness", "HighSB", "HighMeanMAF", "HighMQ0", "InVarRegion", "NonSegregating", "ExcessiveNoCalls"))
#  finalSNPnumbersMatrix <- recombinationPlotSeries(vcfFinalFiltered, plotFilestem=paste(cross, "final", sep="."), filters=c("ExcessiveNoCalls", "HighMissingness", "HighMeanMAF", "HighSB", "HighMQ0", "InVarRegion", "MultiAllelic", "NonSegregating", "HighMaxParentMAF", "HighMaxMAF"), sampleIDcolumn=sampleIDcolumn, sampleIDmappingsColumn=sampleIDmappingsColumn, sampleDuplicates=finalSampleQCresults[["sampleDuplicates"]])
##  qcPlusUniqueSamples <- setdiff(initialSampleQCresults[["uniqueSamples"]], initialSampleQCresults[["qcFailedSamples"]])
##  qcPlusUniqueSamples <- setdiff(finalSampleQCresults[["uniqueSamples"]], finalSampleQCresults[["qcFailedSamples"]])
#  qcPlusUniqueSamples <- setdiff(initialSampleQCresults[["uniqueSamples"]], c(finalSampleQCresults[["qcFailedSamples"]], initialSampleQCresults[["qcFailedSamples"]]))
#  if(file.exists(vcfUniqueFilteredRda) & shouldUseSavedVersions) {
#    load(vcfUniqueFilteredRda)
#  } else {
#    vcfUniqueSamples <- annotateVcf(filterVcf(vcfVariant[, qcPlusUniqueSamples]))
#    qcFilteringResultsUnique <- qcFilteringPlots(vcfUniqueSamples, plotFilestem=paste(cross, "unqiuePreFiltering", sep="."))
#    filt(vcfUniqueSamples) <- "PASS"
#    vcfUniqueFiltered <- setVcfFilters(
#      vcfUniqueSamples,
##      additionalInfoFilters = NULL,
#      additionalInfoFilters = list(
##        "LowQD" = list(column="QD", operator="<=", value=1),
#        "HighSB" = list(column="SB", operator=">=", value=-4200),
#        "HighMeanMAF" = list(column="meanMAF", operator=">=", value=0.1, filterOutNAs=TRUE),
#        "HighMaxMAF" = list(column="maxMAF", operator=">=", value=0.1, filterOutNAs=TRUE),
#        "HighMaxParentMAF" = list(column="maxParentMAF", operator=">=", value=0.1, filterOutNAs=TRUE),
#        "HighMissingness" = list(column="missingness2", operator=">=", value=1, filterOutNAs=TRUE),
#        "HighMQ0" = list(column="MQ0", operator=">=", value=10000, filterOutNAs=TRUE)
#      ),
#      regionsMask                 = varRegions_v3(),
#      shouldSetMultiallelicFilter = TRUE,
#      shouldSetNonSegregatingFilt = TRUE,
#      shouldSetMaxNoCallsFilter   = TRUE
#    )
#    save(vcfUniqueFiltered, file=vcfUniqueFilteredRda)
#  }
#  qcFilteringResultsUniquePostFiltering <- qcFilteringPlots(
#    filterVcf(
#      vcfUniqueFiltered,
#      shouldRemoveInvariant=FALSE,
#      filtersToRemove=c("LowQD", "HighSB", "HighMeanMAF", "HighMissingness", "HighMQ0", "InVarRegion", "ExcessiveNoCalls", "HighMaxParentMAF", "HighMaxMAF")
#    ),
#    plotFilestem=paste(cross, "uniquePostFiltering", sep=".")
#  )
##  finalSNPnumbersMatrix2 <- recombinationPlotSeries(vcfFinalFiltered[, qcPlusUniqueSamples], plotFilestem=paste(cross, "uniqueSamples", sep="."), filters=c("LowQD", "InVarRegion", "NonSegregating", "ExcessiveNoCalls"))
#  uniqueSNPnumbersMatrix <- recombinationPlotSeries(vcfFinalFiltered[, qcPlusUniqueSamples], plotFilestem=paste(cross, "uniqueSamples", sep="."), filters=c("ExcessiveNoCalls", "HighMissingness", "HighMeanMAF", "HighSB", "HighMQ0", "InVarRegion", "MultiAllelic", "NonSegregating", "HighMaxParentMAF", "HighMaxMAF"), sampleIDcolumn=sampleIDcolumn, sampleIDmappingsColumn=sampleIDmappingsColumn, sampleDuplicates=finalSampleQCresults[["sampleDuplicates"]])
#  finalUniqueSampleQCresults <- sampleQC(vcfFinalFiltered[, qcPlusUniqueSamples], discordanceThreshold=discordanceThresholdInitial, plotFilestem=paste(cross, "uniqueSamples", sep="."), gffGRL=gffGRL)
#  vcfSegregating <- vcfFinalFiltered
  
#    vcf <- vcfVariantAnnotated
#    regionsToMask               = varRegions_v3()
#    vcf <- vcf[!(rowData(vcf) %in% regionsToMask)]
#    vcf <- vcf[elementLengths(alt(vcf)) == 1]
#    quantiles <- cut_number(values(info(vcf))[["QD"]], 100)
#    proportions <- by(values(info(vcf))[["MendelianErrors"]], quantiles, function(x) length(which(x>0))/length(x))
  
#  qcFailedSamples <- determineQCfailedSites(
#    vcf,
#    heterozygosityThreshold=heterozygosityThreshold,
#    missingnessThreshold=heterozygosityThreshold
#  )
#  vcfVariantQCplus <- qcPlusVCF(vcfVariant, qcFailedSamples=qcFailedSamples)
#  vcfVariantQCplus <- vcfVariant
#  rm(vcfVariant)
  gc()
#  if(is.null(plotFilestem)) {
#    samplesToUse <- uniqueSamples(vcfFiltered, discordanceThreshold=discordanceThresholdMg) # should output heatmap of discordances
#  } else {
#    samplesToUse <- uniqueSamples(vcfFiltered, discordanceThreshold=discordanceThresholdMg, plotFilestem=plotFilestem) # should output heatmap of discordances
#  }
#  gc()
#  vcfListSegregating <- segregatingSitesList(vcfList, samplesToUse=samplesToUse, pdfFilestem=plotFilestem) #allPaintingSeries, also will probably extend to remove short haplotypes with recombination points in many samples
#  vcfSegregating <- combineVcfListIntoVcf(vcfListSegregating)
#  uniqueSamples(vcfSegregating, discordanceThreshold=discordanceThresholdSeg, plotFilestem=paste(meta(exptData(vcf)[["header"]])["DataSetName", "Value"], "segregating", sep="."))
  if(shouldCompareWithJiang) {
    jiangVcf <- loadJiangGenotypesAsVcf() # make this return a VCF?
    jiangSampleQCresults <- sampleQC(jiangVcf, discordanceThresholdJiang, plotFilestem="JiangEtAl", GTsToIntMapping = c("7"=1, "G"=2, "."=0), shouldCalcMissingnessAndHet=FALSE, shouldRenameSamples=FALSE) # should output heatmap of discordances
#    jiangUniqueSamples <- uniqueSamples(jiangVcf, discordanceThresholdJiang, plotFilestem="JiangEtAl", GTsToIntMapping = c("7"=1, "G"=2, "."=0)) # should output heatmap of discordances
    if(cross=="v3UG_7g8xGb4") {
      seqlevels(jiangVcf) <- sprintf("Pf3D7_%02d_v3", as.integer(sub("MAL", "", seqlevels(jiangVcf))))
    }
    genotypeConcordanceRaw <- compareCalls(vcfVariantAnnotated, jiangVcf, plotFilestem=paste(cross, "comparison", "raw", sep="."), discordanceThreshold=discordanceThresholdRawVsJia)
#    genotypeConcordanceRaw <- compareCalls(vcfFiltered, jiangVcf, plotFilestem=paste(cross, "comparison", "raw", sep="."), discordanceThreshold=discordanceThresholdRawVsJia)
    genotypeConcordance <- compareCalls(vcfSegregating, jiangVcf, plotFilestem=paste(cross, "comparison", "filtered", sep="."), discordanceThreshold=discordanceThresholdFltVsJia) # Should give slide 3, histogram of pair-wise numbers of discordant, heatmap of sample discordances and heatmap for discordances for presumed identical, recombinationPlot of both together
#    genotypeConcordancePf3D7_02_v3 <- compareCalls(vcfSegregating[seqnames(vcfSegregating)=="Pf3D7_02_v3"], jiangVcf[seqnames(jiangVcf)=="Pf3D7_02_v3"], plotFilestem=paste(cross, "comparison", "Pf3D7_02_v3", sep="."), discordanceThreshold=discordanceThreshold5ChVsJia) # Should give slide 3, histogram of pair-wise numbers of discordant, heatmap of sample discordances and heatmap for discordances for presumed identical, recombinationPlot of both together
#    genotypeConcordance6chromosomes <- compareCalls(vcfSegregating[seqnames(vcfSegregating) %in% c("Pf3D7_02_v3", "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3")], jiangVcf[seqnames(jiangVcf) %in% c("Pf3D7_02_v3", "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3")], plotFilestem=paste(cross, "comparison", "6chromosomes", sep="."), discordanceThreshold=discordanceThreshold5ChVsJia) # Should give slide 3, histogram of pair-wise numbers of discordant, heatmap of sample discordances and heatmap for discordances for presumed identical, recombinationPlot of both together
#    genotypeConcordance5chromosomes <- compareCalls(vcfSegregating[seqnames(vcfSegregating) %in% c("Pf3D7_02_v3", "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_10_v3", "Pf3D7_11_v3")], jiangVcf[seqnames(jiangVcf) %in% c("Pf3D7_02_v3", "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_10_v3", "Pf3D7_11_v3")], plotFilestem=paste(cross, "comparison", "5chromosomes", sep="."), discordanceThreshold=discordanceThreshold5ChVsJia) # Should give slide 3, histogram of pair-wise numbers of discordant, heatmap of sample discordances and heatmap for discordances for presumed identical, recombinationPlot of both together
    gc()
    if(!file.exists(paste(cross, "mgRecombinations.rda", sep=".")) | !shouldUseSavedVersions) {
      mgRecombinations <- recombinationPoints(vcfSegregating[filt(vcfSegregating)=="PASS", qcPlusUniqueSamples], gffGRL) # extend crossoversAnalysis to include classification as exonic, intronic, etc
      save(mgRecombinations, file=paste(cross, "mgRecombinations.rda", sep="."))
    } else {
      load(paste(cross, "mgRecombinations.rda", sep="."))
    }
    if(!file.exists("~/jiangRecombinations.rda") | !shouldUseSavedVersions) {
      jiangRecombinations <- recombinationPoints(jiangVcf, gffGRL, GTsToIntMapping = c("7"=1, "G"=2, "."=0))
      save(jiangRecombinations, file="jiangRecombinations.rda")
    } else {
      load("~/jiangRecombinations.rda")
    }
    medianBreakpointAccuracies <- compareRecombinations(mgRecombinations, jiangRecombinations) # to include slide 9 plot, slide 12 plot, median accuracies, Venn
  } else {
    if(!file.exists(paste(cross, "mgRecombinations.rda", sep=".")) | !shouldUseSavedVersions) {
      mgRecombinations <- recombinationPoints(vcfSegregating[filt(vcfSegregating)=="PASS"], gffGRL) # extend crossoversAnalysis to include classification as exonic, intronic, etc
      save(mgRecombinations, file=paste(cross, "mgRecombinations.rda", sep="."))
    } else {
      load(paste(cross, "mgRecombinations.rda", sep="."))
    }
  }
  recombinationRates <- analyseRecombinations(mgRecombinations, plotFilestem=cross) # to include slide 13 plot, breakdown of CO and GC by progeny, chromosome and by cross, CO and GC rates
  returnList <- list(
    vcfVariantAnnotated                  = vcfVariantAnnotated,
    vcfSegregating                       = vcfSegregating,
    mgRecombinations                     = mgRecombinations,
    recombinationRates                   = recombinationRates,
    qcFilteringResults_raw               = qcFilteringResults_raw,
    qcFilteringResults_core              = qcFilteringResults_core,
#    qcFilteringResults                   = qcFilteringResults,
    initialSampleQCresults               = initialSampleQCresults,
    initialSNPnumbersMatrix              = initialSNPnumbersMatrix,
#    qcFilteringResultsFinalPostFiltering = qcFilteringResultsFinalPostFiltering,
#    finalSampleQCresults                 = finalSampleQCresults,
#    finalSNPnumbersMatrix                = finalSNPnumbersMatrix,
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
  save(returnList, file=paste(cross, "returnList.rda", sep="."))
  return(returnList)
}

#Post-hoc debugging
#temp_finalSNPnumbersMatrix <- recombinationPlotSeries(
#  v3UG_3d7xHb3[["vcfSegregating"]],
#  plotFilestem=paste("temp", sep="."),
#  filters=c("ExcessiveNoCalls", "HighMissingness", "HighMeanMAF", "HighSB", "HighMQ0", "InVarRegion", "MultiAllelic", "NonSegregating"),
#  sampleIDcolumn="ena_run_accession",
#  sampleIDmappingsColumn="ena_run_accession",
#  sampleDuplicates=v3UG_3d7xHb3[["finalSampleQCresults"]][["sampleDuplicates"]]
#)

