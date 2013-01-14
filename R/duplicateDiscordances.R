# duplicateDiscordances.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


duplicateDiscordances <- function(
  vcf,
  duplicatesDF,
  possibleMissingValues       = c(".", "./.", ".|.")
) {
  pairsToUse <- which(duplicatesDF[["sampleID1"]] %in% dimnames(vcf)[[2]] & duplicatesDF[["sampleID2"]] %in% dimnames(vcf)[[2]])
  if(length(pairsToUse) == 0) {
    return(0)
  }
  GTsFirstSample <- geno(vcf)[["GT"]][, duplicatesDF[pairsToUse, "sampleID1"]]
  GTsSecondSample <- geno(vcf)[["GT"]][, duplicatesDF[pairsToUse, "sampleID2"]]
  GTsFirstSample[GTsFirstSample %in% possibleMissingValues] <- NA
  GTsSecondSample[GTsSecondSample %in% possibleMissingValues] <- NA
  discordanceMatrix <- GTsFirstSample != GTsSecondSample
}
