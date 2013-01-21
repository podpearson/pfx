# readAllChromosomesWithRegionsMask.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


readAllChromosomesWithRegionsMask <- function(
  cross                       = "3d7_hb3",
  variantType                 = "snps",
  genotypesDirectory          = "data/3d7_v3/bwa_n0.01_k4_l32/genotypes/gatk",
  chromosomes                 = sprintf("Pf3D7_%02d_v3", 1:14),
  genoToLoad                  = c("GT", "AD", "GQ", "MQ0"),
  nonVarRegions               = nonVarRegions_v3(),
  outputDirectory             = sub("genotypes", "genotypes_analysis_20130121", genotypesDirectory),
  genotypesFileFmt            = "%s.annotated.vcf",
  overwriteExisting           = NULL,
  parentalStrains             = NULL
)  {
  if(!file.exists(outputDirectory)) {
    dir.create(outputDirectory, recursive=TRUE)
  }
  if(is.null(overwriteExisting)) {
    lapply(
      chromosomes,
      function(chromosome) {
        inputVcfFilename <- file.path(genotypesDirectory, cross, variantType, sprintf(genotypesFileFmt, chromosome))
        outputVcfFilename <- file.path(outputDirectory, cross, variantType, sprintf(genotypesFileFmt, chromosome))
        outputVcfGzFilename <- paste(outputVcfFilename, "gz", sep=".")
        outputRdaFilename <- paste(outputVcfFilename, "coreAllSamples.rda", sep=".")
        if(!file.exists(outputVcfGzFilename)) {
          if(!file.exists(outputVcfFilename)) {
            file.copy(inputVcfFilename, outputVcfFilename)
          }
          bgzip(outputVcfFilename)
          indexTabix(outputVcfGzFilename, "vcf4")
        }
        readSingleChromosomeWithRegionsMask(
          vcfFilename                 = outputVcfGzFilename,
          chromosome                  = chromosome,
          genoToLoad                  = genoToLoad,
          nonVarRegions               = nonVarRegions,
          parentalStrains             = parentalStrains,
          outputRdaFilename           = outputRdaFilename
        )
      }
    )
  } else {
    lapply(
      chromosomes,
      function(chromosome) {
        readSingleChromosomeWithRegionsMask(
          file.path(genotypesDirectory, cross, variantType, sprintf(genotypesFileFmt, chromosome)),
          chromosome,
          filtersToRemove             = filtersToRemove,
          samplesToRemove             = samplesToRemove,
          overwriteExisting           = overwriteExisting,
          parentalStrains             = parentalStrains,
          outputRdaFilename           = file.path(outputDirectory, cross, variantType, sub("\\.vcf[\\.gz]*$", paste("\\.variantSites\\.", chromosome,"\\.rda", sep=""), sprintf(genotypesFileFmt, chromosome)))
        )
      }
    )
  }
}

