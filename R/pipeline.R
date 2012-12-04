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

#v3UG_3d7xHb3 <- pipeline("v3UG_3d7xHb3", "/data/malariagen2/plasmodium/pf-crosses/data/3d7_v3/bwa_default/snp_genotypes_analysis/PFproj2-unigen-snponly-BQ20.vcf.gz", chromosomes=sprintf("Pf3D7_%02d_v3", 1:14), filtersToRemove="LowQual", overwriteExisting=TRUE, discordanceThresholdMg=28000, plotFilestem="v3UG_3d7xHb3")

pipeline <- function(
  cross                       = "7g8xGb4",
  vcfFilename                 = file.path("/data/galton/users/rpearson", "crossesTesting", "release", paste(cross, "-qcPlusSamples-0.1.vcf.gz", sep="")),
  gffFilename                 = "/data/galton/mirror/nfs/team112/annotation/plasmodium/falciparum/Pfalciparum_PlasmoDB-7.2.gff",
  gffGRL                      = readGffAsGRangesList(gffFilename),
  chromosomes                 = sprintf("MAL%d", 1:14),
#  chromosomes                 = sprintf("Pf3D7_%02d_v3", 1:14),
  filtersToRemove             = NULL,
  samplesToRemove             = NULL,
#  samplesToRemove             = c("ERR045643", "ERR045644", "ERR045645", "ERR045646", "ERR045647", "ERR045625"),
  plotFilestem                = NULL,
  heterozygosityThreshold     = 0.3,
  missingnessThreshold        = 0.3,
  discordanceThresholdMg      = 10000,
  discordanceThresholdSeg     = 100,
  discordanceThresholdJiang   = 100,
  discordanceThresholdRawVsJia= 100,
  discordanceThresholdFltVsJia= 100,
  shouldCompareWithJiang      = grepl("7g8xGb4", cross),
  overwriteExisting           = NULL,
  shouldUseSavedVersions      = TRUE,
  vcfListRda                  = paste(cross, "vcfList.rda", sep="."),
  vcfVariantRda               = paste(cross, "vcfVariant.rda", sep="."),
  vcfVariantAnnotatedRda      = paste(cross, "vcfVariantAnnotated.rda", sep=".")
) {
  if(file.exists(vcfListRda) & shouldUseSavedVersions) {
    load(vcfListRda)
  } else {
    if(is.null(overwriteExisting)) {
      vcfList <- variantSitesRdaFiles(vcfFilename, chromosomes=chromosomes, filtersToRemove=filtersToRemove, samplesToRemove=samplesToRemove)
    } else {
      vcfList <- variantSitesRdaFiles(vcfFilename, chromosomes=chromosomes, filtersToRemove=filtersToRemove, samplesToRemove=samplesToRemove, overwriteExisting=overwriteExisting)
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
    vcfVariantAnnotated <- annotateSegregationStatus(vcfVariant)
    save(vcfVariantAnnotated, file=vcfVariantAnnotatedRda)
  }
  qcFilteringResults <- qcFilteringPlots(vcfVariantAnnotated, plotFilestem=cross)
#  browser()
  
#  qcFailedSamples <- determineQCfailedSites(
#    vcf,
#    heterozygosityThreshold=heterozygosityThreshold,
#    missingnessThreshold=heterozygosityThreshold
#  )
#  vcfVariantQCplus <- qcPlusVCF(vcfVariant, qcFailedSamples=qcFailedSamples)
  vcfVariantQCplus <- vcfVariant
  rm(vcfVariant)
  gc()
  if(is.null(plotFilestem)) {
    samplesToUse <- uniqueSamples(vcfVariantQCplus, discordanceThreshold=discordanceThresholdMg) # should output heatmap of discordances
  } else {
    samplesToUse <- uniqueSamples(vcfVariantQCplus, discordanceThreshold=discordanceThresholdMg, plotFilestem=plotFilestem) # should output heatmap of discordances
  }
  gc()
  vcfListSegregating <- segregatingSitesList(vcfList, samplesToUse=samplesToUse, pdfFilestem=plotFilestem) #allPaintingSeries, also will probably extend to remove short haplotypes with recombination points in many samples
  vcfSegregating <- combineVcfListIntoVcf(vcfListSegregating)
#  uniqueSamples(vcfSegregating, discordanceThreshold=discordanceThresholdSeg, plotFilestem=paste(meta(exptData(vcf)[["header"]])["DataSetName", "Value"], "segregating", sep="."))
  if(shouldCompareWithJiang) {
    jiangVcf <- loadJiangGenotypesAsVcf() # make this return a VCF?
    jiangUniqueSamples <- uniqueSamples(jiangVcf, discordanceThresholdJiang, plotFilestem="JiangEtAl", GTsToIntMapping = c("7"=1, "G"=2, "."=0)) # should output heatmap of discordances
    if(cross=="v3UG_7g8xGb4") {
      seqlevels(jiangVcf) <- sprintf("Pf3D7_%02d_v3", as.integer(sub("MAL", "", seqlevels(jiangVcf))))
    }
    genotypeConcordanceRaw <- compareCalls(vcfVariantQCplus, jiangVcf, plotFilestem=paste(cross, "comparison", "raw", sep="."), discordanceThreshold=discordanceThresholdRawVsJia)
    genotypeConcordance <- compareCalls(vcfSegregating, jiangVcf, plotFilestem=paste(cross, "comparison", "filtered", sep="."), discordanceThreshold=discordanceThresholdFltVsJia) # Should give slide 3, histogram of pair-wise numbers of discordant, heatmap of sample discordances and heatmap for discordances for presumed identical, recombinationPlot of both together
    gc()
    if(!file.exists(paste(cross, "mgRecombinations.rda", sep="."))) {
      mgRecombinations <- recombinationPoints(vcfSegregating, gffGRL) # extend crossoversAnalysis to include classification as exonic, intronic, etc
      save(mgRecombinations, file=paste(cross, "mgRecombinations.rda", sep="."))
    } else {
      load(paste(cross, "mgRecombinations.rda", sep="."))
    }
    if(!file.exists("~/jiangRecombinations.rda")) {
      jiangRecombinations <- recombinationPoints(jiangVcf, gffGRL, GTsToIntMapping = c("7"=1, "G"=2, "."=0))
      save(jiangRecombinations, file="jiangRecombinations.rda")
    } else {
      load("~/jiangRecombinations.rda")
    }
    medianBreakpointAccuracies <- compareRecombinations(mgRecombinations, jiangRecombinations) # to include slide 9 plot, slide 12 plot, median accuracies, Venn
  } else {
    if(!file.exists(paste(cross, "mgRecombinations.rda", sep="."))) {
      mgRecombinations <- recombinationPoints(vcfSegregating, gffGRL) # extend crossoversAnalysis to include classification as exonic, intronic, etc
      save(mgRecombinations, file=paste(cross, "mgRecombinations.rda", sep="."))
    } else {
      load(paste(cross, "mgRecombinations.rda", sep="."))
    }
  }
  recombinationRates <- analyseRecombinations(mgRecombinations, plotFilestem=cross) # to include slide 13 plot, breakdown of CO and GC by progeny, chromosome and by cross, CO and GC rates
  if(shouldCompareWithJiang) {
    return(list(vcfSegregating=vcfSegregating, mgRecombinations=mgRecombinations, jiangRecombinations=jiangRecombinations, genotypeConcordance=genotypeConcordance, medianBreakpointAccuracies=medianBreakpointAccuracies, recombinationRates=recombinationRates))
  } else {
    return(list(vcfSegregating=vcfSegregating, mgRecombinations=mgRecombinations, recombinationRates=recombinationRates))
  }
}
