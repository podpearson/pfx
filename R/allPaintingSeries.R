# allPaintingSeries.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################

#temp_7g8xGb4 <- allPaintingSeries(file.path("/data/galton/users/rpearson", "crossesTesting", "release", paste("7g8xGb4-qcPlusSamples-0.1.vcf.MAL", 1:14, ".rda", sep="")))
#temp_3d7xHb3 <- allPaintingSeries(file.path("/data/galton/users/rpearson", "crossesTesting", "release", paste("3d7xHb3-qcPlusSamples-0.1.vcf.MAL", 1:14, ".rda", sep="")))
#temp_Hb3xDd2 <- allPaintingSeries(file.path("/data/galton/users/rpearson", "crossesTesting", "release", paste("Hb3xDd2-qcPlusSamples-0.1.vcf.MAL", 1:14, ".rda", sep="")))

allPaintingSeries <- function(
  vcfRdaFilenames               = file.path("/data/galton/users/rpearson", "crossesTesting", "release", paste("7g8xGb4-qcPlusSamples-0.1.vcf.MAL", 1:14, ".rda", sep="")),
  vcfToUseForSampleChecking     = 4,
  filtersToRemove               = c("NoAlternative"),
  discordanceThreshold          = 500,
  verbose                       = TRUE
) {
   if(verbose) {
    cat("allPaintingSeries: loading vcf...")
  }
  load(vcfRdaFilenames[vcfToUseForSampleChecking])
  if(verbose) {
    cat("done\n")
  }
  if(verbose) {
    cat("allPaintingSeries: filtering vcf...")
  }
  sapply(
    filtersToRemove,
    function(filterToRemove) {
      vcf <- vcf[!grepl(filterToRemove, filt(vcf))]
    }
  )
  samplesToUse <- uniqueSamples(vcf, discordanceThreshold=discordanceThreshold)
  lapply(
    vcfRdaFilenames,
    function(vcfRdaFilename) {
      paintingSeries(vcfRdaFilename, samplesToUse=samplesToUse, verbose=verbose)
    }
  )
}
