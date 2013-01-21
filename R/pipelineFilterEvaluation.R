# pipelineFilterEvaluation.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


# pipelineFilterEvaluation(chromosomes = sprintf("Pf3D7_%02d_v3", 1))

pipelineFilterEvaluation <- function(
  cross                       = "3d7_hb3",
  variantType                 = "snps",
  genotypesDirectory          = "data/3d7_v3/bwa_n0.01_k4_l32/genotypes/per_sample_realigned/gatk",
  chromosomes                 = sprintf("Pf3D7_%02d_v3", 1:14),
  outputDirectory             = "data/3d7_v3/bwa_n0.01_k4_l32/genotypes_analysis_20130121/per_sample_realigned/gatk",
#  outputDirectory             = sub("genotypes", "genotypes_analysis", genotypesDirectory),
  genotypesFileFmt            = "%s.annotated.vcf",
  gffFilename                 = "/data/malariagen2/plasmodium/pf-crosses/data/genome/sanger/version3/September_2012/Pf3D7_v3.gatk.gff",
  gffGRL                      = readGffAsGRangesList(gffFilename, chromsomeNames=chromosomes),
  parentalStrains             = if(cross=="3d7_hb3") {
    c("ERR019061", "ERR019054")
  } else if(cross=="7g8_gb4") {
    c("ERR027099", "ERR027100")
  } else if(cross=="hb3_dd2") {
    c("ERR012788", "ERR012840")
  },
  sampleIDcolumn              = "ena_run_accession",
  sampleIDmappingsColumn      = sampleIDcolumn,
  filtersToRemove             = NULL,
  samplesToRemove             = NULL,
#  samplesToRemove             = c("ERR045643", "ERR045644", "ERR045645", "ERR045646", "ERR045647", "ERR045625"),
  plotFilestem                = NULL,
#  heterozygosityThreshold     = 0.3,
#  missingnessThreshold        = 0.3,
#  discordanceThresholdInitial = 1000,
#  discordanceThresholdMg      = 10000,
#  discordanceThresholdSeg     = 100,
#  discordanceThresholdJiang   = 100,
#  discordanceThresholdRawVsJia= 100,
#  discordanceThresholdFltVsJia= 100,
#  discordanceThreshold5ChVsJia= 50,
#  shouldCompareWithJiang      = grepl("7g8xGb4", cross),
  overwriteExisting           = NULL,
  shouldUseSavedVersions      = TRUE,
  vcfCoreAllSamplesRda        = file.path(outputDirectory, cross, variantType, paste(cross, "vcfCoreFinalSamples.rda", sep="."))
#  vcfListRda                  = file.path(outputDirectory, cross, variantType, paste(cross, "vcfList.rda", sep=".")),
#  vcfVariantRda               = file.path(outputDirectory, cross, variantType, paste(cross, "vcfVariant.rda", sep=".")),
#  vcfVariantAnnotatedRda      = file.path(outputDirectory, cross, variantType, paste(cross, "vcfVariantAnnotated.rda", sep=".")),
#  vcfInitialFilteredRda       = file.path(outputDirectory, cross, variantType, paste(cross, "vcfInitialFiltered.rda", sep=".")),
#  initialSampleQCresultsRda   = file.path(outputDirectory, cross, variantType, paste(cross, "initialSampleQCresults.rda", sep=".")),
#  vcfFinalFilteredRda         = file.path(outputDirectory, cross, variantType, paste(cross, "vcfFinalFiltered.rda", sep=".")),
#  vcfUnfilteredFinalSamplesRda = file.path(outputDirectory, cross, variantType, paste(cross, "vcfUnfilteredFinalSamples.rda", sep=".")),
#  vcfCoreFinalSamplesRda      = file.path(outputDirectory, cross, variantType, paste(cross, "vcfCoreFinalSamples.rda", sep=".")),
#  vcfUniqueFilteredRda        = file.path(outputDirectory, cross, variantType, paste(cross, "vcfUniqueFiltered.rda", sep="."))
) {
  if(file.exists(vcfCoreAllSamplesRda) & shouldUseSavedVersions) {
    load(vcfCoreAllSamplesRda)
  } else {
    if(is.null(overwriteExisting)) {
      vcfList <- readAllChromosomesWithRegionsMask(
        cross                       = cross,
        variantType                 = variantType,
        genotypesDirectory          = genotypesDirectory,
        chromosomes                 = chromosomes,
        outputDirectory             = outputDirectory,
        genotypesFileFmt            = genotypesFileFmt,
        parentalStrains             = parentalStrains
      )
    } else {
      vcfList <- readAllChromosomesWithRegionsMask(
        cross                       = cross,
        variantType                 = variantType,
        genotypesDirectory          = genotypesDirectory,
        chromosomes                 = chromosomes,
        outputDirectory             = outputDirectory,
        genotypesFileFmt            = genotypesFileFmt,
        overwriteExisting           = overwriteExisting,
        parentalStrains             = parentalStrains
      )
    }
    vcfCoreAllSamples <- combineVcfListIntoVcf(vcfList)
    save(vcfCoreAllSamples, file=vcfCoreAllSamplesRda)
  }
}

