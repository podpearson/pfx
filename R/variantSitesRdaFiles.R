# variantSitesRdaFiles.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


variantSitesRdaFiles <- function(
  vcfFilename                 = file.path("/data/galton/users/rpearson", "crossesTesting", "release", "7g8xGb4-qcPlusSamples-0.1.vcf.gz"),
  chromosomes                 = c(paste("MAL", 1:14, sep="")),
#  chromosomes                 = c(paste("MAL", 1:14, sep=""), "MITO1", paste("PLAST", 1:2, sep=""))
  filtersToRemove             = NULL,
  samplesToRemove             = NULL,
  overwriteExisting           = NULL,
  parentalStrains             = NULL
)  {
  if(is.null(overwriteExisting)) {
    lapply(
      chromosomes,
      function(x) {
        createSingleChromosomeVariantSitesRdaFile(
          vcfFilename,
          x,
          filtersToRemove=filtersToRemove,
          samplesToRemove=samplesToRemove,
          parentalStrains=parentalStrains
        )
      }
    )
  } else {
    lapply(
      chromosomes,
      function(x) {
        createSingleChromosomeVariantSitesRdaFile(
          vcfFilename,
          x,
          filtersToRemove=filtersToRemove,
          samplesToRemove=samplesToRemove,
          overwriteExisting=overwriteExisting,
          parentalStrains=parentalStrains
        )
      }
    )
  }
}
