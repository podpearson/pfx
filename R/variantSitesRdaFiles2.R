# variantSitesRdaFiles2.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


variantSitesRdaFiles2 <- function(
  cross                       = "3d7_hb3",
  variantType                 = "snps",
  genotypesDirectory          = "data/3d7_v3/bwa_n0.01_k4_l32/genotypes/gatk",
  chromosomes                 = sprintf("Pf3D7_%02d_v3", 1:14),
  outputDirectory             = sub("genotypes", "genotypes_analysis", genotypesDirectory),
  genotypesFileFmt            = "%s.raw.vcf",
#  chromosomes                 = c(paste("MAL", 1:14, sep=""), "MITO1", paste("PLAST", 1:2, sep=""))
  filtersToRemove             = NULL,
  samplesToRemove             = NULL,
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
        createSingleChromosomeVariantSitesRdaFile(
          file.path(genotypesDirectory, cross, variantType, sprintf(genotypesFileFmt, chromosome)),
          chromosome,
          filtersToRemove             = filtersToRemove,
          samplesToRemove             = samplesToRemove,
          parentalStrains             = parentalStrains,
          outputRdaFilename           = file.path(outputDirectory, cross, variantType, sub("\\.vcf[\\.gz]*$", paste("\\.variantSites\\.", chromosome,"\\.rda", sep=""), sprintf(genotypesFileFmt, chromosome)))
        )
      }
    )
  } else {
    lapply(
      chromosomes,
      function(chromosome) {
        createSingleChromosomeVariantSitesRdaFile(
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
